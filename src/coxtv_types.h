#ifndef COXTV_TYPES_H
#define COXTV_TYPES_H


#include <Eigen/Dense>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::VectorXi VectorXi;

typedef Eigen::Map<const Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<const Eigen::MatrixXf> MapMatf;
typedef Eigen::Map<const Eigen::VectorXd> MapVecd;
typedef Eigen::Map<const Eigen::VectorXi> MapVeci;
typedef Eigen::Map<const Eigen::MatrixXi> MapMati;
typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> PermMat;

#endif
