#include <Rcpp.h>
#include "coxtv_types.h"
#include <omp.h>
#include <iostream>
#include <sys/time.h>
#include <cmath>


class Cox_prob_td{
    const int p, N, ptd;
    const int num_events, num_unique_events, vtotal;
    MapMatd X; // Here X is N by p
    MapMatd V; // The time-dependent covariates, Vind[-1] by ptd
    MapVeci Vind; // starting cols of each unique event time
    MapVeci order; // order of all observations
    MapVeci erank; // unique rank of the events (min rank)
    MapVeci eind; //event index in X
    MapVecd te_count; //count of tied_events at each unique event rank

    VectorXd eta;
    VectorXd grad_baseline;
    VectorXd eta_td;
    VectorXd exp_val;
    VectorXd A_cumu;
    PermMat rev_order;
    PermMat order_mat;

    VectorXi non_unique_ranks; // event ranks that are not unique
    int num_non_unique_ranks;

    void get_grad_baseline(){
        grad_baseline.resize(ptd+p);
        grad_baseline.segment(0, ptd) = (V(Vind.segment(0, num_unique_events), Eigen::all).colwise().sum()).transpose();
        // Adjust for tied events
        for(int i = 0; i < num_non_unique_ranks; ++i){
            int ind = non_unique_ranks[i];
            for(int j = 1; j < te_count[ind]; ++j){
                grad_baseline.segment(0, ptd) += V.row(Vind[ind]+j);
            }
        }
        grad_baseline.segment(ptd, p) = (X(eind, Eigen::all).colwise().sum()).transpose();
        grad_baseline *= -1.0;
    }



    public:
    Cox_prob_td(int p,
                int N,
                int ptd,
                int vtotal,
                int num_events,
                int num_unique_events,
                const double *X,
                const double *V,
                const int *Vind,
                const int *order,
                const int *erank,
                const int *eind,
                const double *te_count) : p(p),
                                        N(N),
                                        ptd(ptd),
                                        num_events(num_events),
                                        num_unique_events(num_unique_events),
                                        vtotal(vtotal),
                                        X(X, N, p),
                                        V(V, vtotal, ptd),
                                        Vind(Vind, num_unique_events+1),
                                        order(order, N),
                                        erank(erank, num_unique_events),
                                        eind(eind, num_events),
                                        te_count(te_count, num_unique_events),
                                        rev_order(this->order),
                                        non_unique_ranks(num_unique_events)
        {
            eta.resize(N);
            eta_td.resize(vtotal);
            exp_val.resize(N);
            A_cumu.resize(N);
            num_non_unique_ranks = 0;
            for(int i = 0; i < num_unique_events; ++i){
                if(te_count[i] > 1){
                    non_unique_ranks[num_non_unique_ranks] = i;
                    num_non_unique_ranks++;
                }
            }
            VectorXi tmp = non_unique_ranks(Eigen::seq(0, num_non_unique_ranks - 1));
            non_unique_ranks = tmp;
            order_mat = rev_order.transpose();
            get_grad_baseline();
        }



    double get_gradient(const VectorXd & v, VectorXd & grad, bool get_val=false)
    {
        double cox_val = 0.0;
        A_cumu.setZero();
        grad = grad_baseline;
        eta.noalias() =  (X * v.segment(ptd, p)); // X is short and fat
        if(get_val)
        {
            cox_val -= eta(eind).sum();
        }
        eta = order_mat * eta;

        eta_td.noalias() = V * v.segment(0, ptd); // V is tall and skinny
        // unfortunately it seems the following is necessary
        // and not sure now to parallelize without blowing up memory
        int rank, len, td_start;
        double denom;
        for(int i = 0; i < num_unique_events; ++i)
        {
            rank = erank[i];
            len = N - rank; // size of riskset
            td_start = Vind[i];

            // Need to make sure that order and V are aligned properly
            exp_val.segment(0, len) = (eta_td.segment(td_start, len) +  eta.segment(rank, len)).array().exp().matrix();
            denom = exp_val.segment(0, len).sum();
            cox_val += log(denom)*(te_count[i]);
            exp_val.segment(0, len) *= ((te_count[i])/denom);
            A_cumu.segment(rank, len) += exp_val.segment(0, len);
            grad.segment(0, ptd) += (exp_val.segment(0, len).transpose() * V.block(td_start, 0, len, ptd)).transpose();
        }
        grad.segment(ptd, p) += ((rev_order * A_cumu).transpose() * X).transpose();
        grad /= (double)N;

        if(get_val){
            cox_val -= eta_td(Vind.segment(0, num_unique_events)).sum();

            // Adjust for tied events
            for(int i = 0; i < num_non_unique_ranks; ++i){
                int ind = non_unique_ranks[i];
                // for(int j = 1; j < te_count[ind]; ++j){
                //     //cox_val -= eta_td[Vind[ind]+j];
                // }
                cox_val -= eta_td.segment(Vind[ind]+1, te_count[ind]-1).sum();
            }
            cox_val /= (double)N;
        }

        return cox_val;
    }

    double get_value(const VectorXd & v)
    {
        eta.noalias() =  (X * v.segment(ptd, p)); // X is short and fat
        double cox_val = -eta(eind).sum();
        eta = order_mat * eta;

        eta_td.noalias() = V * v.segment(0, ptd); // V is tall and skinny

        // This is parallelizable
        int rank, len, td_start;
        double denom;
        for(int i = 0; i < num_unique_events; ++i)
        {
            rank = erank[i];
            len = N - rank; // size of riskset
            td_start = Vind[i];

            // Need to make sure that order and V are aligned properly
            exp_val.segment(0, len) = (eta_td.segment(td_start, len) +  eta.segment(rank, len)).array().exp().matrix();
            denom = exp_val.segment(0, len).sum();
            cox_val += log(denom)*te_count[i];
        }

        cox_val -= eta_td(Vind.segment(0, num_unique_events)).sum();

        // Adjust for tied events
        for(int i = 0; i < num_non_unique_ranks; ++i){
            int ind = non_unique_ranks[i];
            cox_val -= eta_td.segment(Vind[ind]+1, te_count[ind]-1).sum();
        }
        cox_val /= (double)N;

        return cox_val;
    }

    void get_residual(const MapVecd & v, VectorXd & result){
        result.setZero();
        A_cumu.setZero();
        result(eind) = (result(eind).array() - 1.0).matrix();
        eta.noalias() =  order_mat * (X * v.segment(ptd, p));
        eta_td.noalias() = V * v.segment(0, ptd);

        int rank, len, td_start;
        double denom;
        for(int i = 0; i < num_unique_events; ++i)
        {
            rank = erank[i];
            len = N - rank; // size of riskset
            td_start = Vind[i];

            // Need to make sure that order and V are aligned properly
            exp_val.segment(0, len) = (eta_td.segment(td_start, len) +  eta.segment(rank, len)).array().exp().matrix();
            denom = exp_val.segment(0, len).sum();
            exp_val.segment(0, len) *= ((te_count[i])/denom);
            A_cumu.segment(rank, len) += exp_val.segment(0, len);
        }
        result = result + (rev_order * A_cumu);
        result /= (double)N;
    }
};

void update_parameters(VectorXd & B, const VectorXd & grad, const VectorXd &v, const double step_size,
                       double lambda_1, const MapVecd & penalty_factor)
{
    B.noalias() = v - step_size*grad;
    // Apply proximal operator here:
    //Soft-thresholding
    B = ((B.cwiseAbs() - lambda_1*step_size*penalty_factor).array().max(0) * B.array().sign()).matrix();
}



// [[Rcpp::export]]
Rcpp::List fit(Rcpp::NumericMatrix X,
                Rcpp::NumericMatrix V,
                Rcpp::IntegerVector Vind,
                Rcpp::IntegerVector order,
                Rcpp::IntegerVector erank,
                Rcpp::IntegerVector eind,
                Rcpp::NumericVector te_count,
                double step_size,
                Rcpp::NumericVector B_init, // initial parameter
                Rcpp::NumericVector lambda_1_all,
                Rcpp::NumericVector penalty_factor,
                int niter,
                double linesearch_beta,
                double eps, // convergence criteria,
                double tol=1e-10 // line search tolerance
                )
{

    int N = X.rows();
    int p = X.cols();
    int ptd = V.cols();
    int vtotal = Vind[Vind.size() - 1];
    int num_events = eind.size();
    int num_unique_events = erank.size();
    const int nlambda = lambda_1_all.size();
    MapVecd pfac(&penalty_factor(0), p+ptd);
    MapVecd Bmap(&B_init(0), p+ptd);


    Cox_prob_td prob(p,
                    N,
                    ptd,
                    vtotal,
                    num_events,
                    num_unique_events,
                    &X(0,0),
                    &V(0,0),
                    &Vind(0),
                    &order(0),
                    &erank(0),
                    &eind(0),
                    &te_count(0));

    VectorXd B(Bmap);
    VectorXd prev_B(B);
    VectorXd v(B);
    VectorXd grad(p+ptd);
    VectorXd grad_ls(p+ptd);

    double cox_val;
    double cox_val_next;
    double rhs_ls; // right-hand side of line search condition
    double lambda_1;
    double step_size_intial = step_size;
    Rcpp::List result(nlambda);
    bool stop; // Stop line searching
    double weight_old, weight_new;
    struct timeval start, end;

    for (int lam_ind = 0; lam_ind < nlambda; ++lam_ind){
        gettimeofday(&start, NULL);

        lambda_1 = lambda_1_all[lam_ind];
        weight_old = 1.0;
        step_size = step_size_intial;

        for (int i = 0; i< niter; i++){
            prev_B.noalias() = B;
            cox_val = prob.get_gradient(v, grad, true);
            while (true){
                update_parameters(B, grad, v, step_size, lambda_1, pfac);

                cox_val_next = prob.get_value(B);
                stop = false;
                if(abs((cox_val_next - cox_val)/fmax(1.0, abs(cox_val_next))) > tol){
                    rhs_ls = cox_val + (grad.array() * (B - v).array()).sum() + (B-v).squaredNorm()/(2*step_size);
                    stop = (cox_val_next <= rhs_ls);
                } else {
                    prob.get_gradient(B, grad_ls, false);
                    rhs_ls = ((B-v).array() * (grad_ls - grad).array()).sum();
                    stop = (abs(rhs_ls) <= (B-v).squaredNorm()/(2*step_size));
                }

                if (stop){
                    break;
                }
                step_size /= linesearch_beta;
                std::cout << step_size << std::endl;
                if(step_size < 1e-5){
                    return result;
                }

            }

            if((prev_B - B).lpNorm<Eigen::Infinity>() < eps){
                std::cout << "convergence based on parameter change reached in " << i <<" iterations\n";
                std::cout << "current step size is " << step_size << std::endl;
                gettimeofday(&end, NULL);
                double delta  = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                std::cout <<  "elapsed time is " << delta << " seconds" << std::endl;
                Rcpp::checkUserInterrupt();
                break;
            }

            // Nesterov weight
            weight_new = 0.5*(1+sqrt(1+4*weight_old*weight_old));
            v.noalias() = B + ((weight_old - 1)/weight_new) * (B - prev_B);
            weight_old = weight_new;
            // v.noalias() = B + ((double)i/(double)(i+3))*(B - previous_B); // Another strategy
            if (i != 0 && i % 100 == 0){
                std::cout << "reached " << i << " iterations\n";
                gettimeofday(&end, NULL);
                double delta  = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                std::cout <<  "elapsed time is " << delta  << " seconds" << std::endl;
                Rcpp::checkUserInterrupt();
            }

        }
        Rcpp::NumericVector RB(B.begin(), B.end());
        result[lam_ind] = RB;
        std::cout << "Solution for the " <<  lam_ind+1 << "th lambda pair is obtained\n";
    }
    return result;
}




// [[Rcpp::export]]
Rcpp::IntegerVector get_te_count(Rcpp::IntegerVector erank){
    Rcpp::IntegerVector result(erank.size());
    int last = -1;
    int ind = -1;
    for(int i = 0; i < erank.size(); ++i){
        if(erank[i] == last){
            result[ind]++;
        } else {
            last = erank[i];
            ind++;
            result[ind] = 1;
        }
    }

    return result[Rcpp::Range(0, ind)];

}


// [[Rcpp::export]]
Rcpp::NumericVector get_residual(Rcpp::NumericMatrix X,
                                Rcpp::NumericMatrix V,
                                Rcpp::IntegerVector Vind,
                                Rcpp::IntegerVector order,
                                Rcpp::IntegerVector erank,
                                Rcpp::IntegerVector eind,
                                Rcpp::NumericVector te_count,
                                Rcpp::NumericVector beta)
{

    int N = X.rows();
    int p = X.cols();
    int ptd = V.cols();
    int vtotal = Vind[Vind.size() - 1];
    int num_events = eind.size();
    int num_unique_events = erank.size();
    MapVecd betamap(&beta(0), p+ptd);

    Cox_prob_td prob(p,
                    N,
                    ptd,
                    vtotal,
                    num_events,
                    num_unique_events,
                    &X(0,0),
                    &V(0,0),
                    &Vind(0),
                    &order(0),
                    &erank(0),
                    &eind(0),
                    &te_count(0));
    VectorXd residual(N);
    prob.get_residual(betamap, residual);
    Rcpp::NumericVector result(residual.begin(), residual.end());
    return result;
}

// Do some of the work in R
// [[Rcpp::export]]
double CCIndex(Rcpp::NumericVector eta, // eta in the right order
                Rcpp::NumericVector eta_td, // Same order as before, indexing using vind
                Rcpp::IntegerVector Vind,
                Rcpp::IntegerVector rank_all, // rank of everyone
                Rcpp::IntegerVector erankwt, // erankwt can have ties
                Rcpp::IntegerVector eranke, // events' rank among events, no increment at ties
                Rcpp::IntegerVector te_count
)
{
    int num_events = erankwt.size();
    int N = eta.size();
    int C = 0;
    int D = 0;
    int T = 0;
    int P = 0; // tied prediction

    #pragma omp parallel for schedule(dynamic, 1) reduction(+:C,D,P,T)
    for(int i = 0; i < num_events;++i){
        int rank = erankwt[i];
        int unique_eind = eranke[i];
        int vstart = Vind[unique_eind];
        int te = te_count[unique_eind]; // Number of tied events at this event time
        double escore = eta[rank] + eta_td[vstart];
        T += (te - 1);
        for(int j = rank+te; j < N; ++j){
            double score = eta[j] + eta_td[vstart+j-rank];
            C += (int)(score < escore);
            D += (int)(score > escore);
            P += (int)(score == escore);
        }
    }
    double denom = (double)C + (double)P + (double)D + ((double)T)/2.0;
    double result = (double)C + ((double)P)/(2.0) + ((double)T)/4.0;
    result /= denom;
    return result;
}


// use 0-based index for every object here
// [[Rcpp::export]]
void get_V(Rcpp::NumericVector V,
           const Rcpp::NumericVector td_covs,
           const Rcpp::NumericVector td_time,
           const Rcpp::NumericVector etime,
           const Rcpp::IntegerVector cumu_measurement,
           const Rcpp::IntegerVector erank,
           const Rcpp::IntegerVector vind)
{
    int ecount = erank.size(); // number of events
    for(int i = 0; i < ecount; ++i){
        int vstart = vind[i];
        int vend = vind[i+1];
        int rank = erank[i];
        double time = etime[i];
        for (int j = vstart; j < vend; ++j){
            int jstart = cumu_measurement[rank];
            int jend = cumu_measurement[rank+1];
            // binary search to find the data in td_covs
            // Find the largest time in td_covs for the jth person that's less than time
            int ind;
            do {
                ind = (jstart + jend)/2;
                if(td_time[ind] >= time){
                    jend = ind;
                } else {
                    jstart = ind;
                }
            } while(jstart < jend - 1);
            V[j] = td_covs[jstart];
            rank++;
        }
    }
}
