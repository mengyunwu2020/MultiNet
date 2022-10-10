#include "math.h"
#include <vector>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using Rcpp::as;
using Eigen::Map;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]
Eigen::MatrixXd Test(NumericMatrix AA);
double thresholdl1(double x, double thr);
//[[Rcpp::export]]
Eigen::ArrayXd sqrt_regression(Eigen::ArrayXd ww, Eigen::MatrixXd X1, Eigen::ArrayXd Y, NumericVector groupLen, NumericVector rangeGroupInd, NumericVector lambda2, double lambda1, int nlambda, int d, int kG, int num_relaxation_round, int max_iter, int innerIter, int standard)
{
  Eigen::ArrayXd Xb, r, w1,gr;
  Eigen::ArrayXXd X=X1;
  int n = X.rows();

  w1.resize(d);
  w1=ww;

  Xb.resize(n);
  Xb=X1*w1.matrix();

  double prec = 1e-4;


  NumericVector meansX;
  meansX=X1.colwise().sum()/n;
  NumericVector varX;
  varX=X1.colwise().sum()/n;
  //
  if(standard==1){

    for(int ii=0;ii<d;ii++){


      varX[ii]=(X.col(ii)- meansX[ii]).square().sum();
      for(int nn=0;nn<n;nn++){
        X1(nn,ii)=(X1(nn,ii)-meansX[ii])/varX[ii];

      }

      w1[ii]=w1[ii]*varX[ii];
    }
    Eigen::ArrayXXd X=X1;
  }



  gr.resize(kG);
  gr.setZero();




  r.resize(n);
  r.setZero();

  Eigen::ArrayXd Xb_master(n);
  Eigen::ArrayXd w1_master(d);


  std::vector<int> actset_group(kG, 0);
  std::vector<int> actset_group_master(kG, 0);



  std::vector<int> actset_idx;
  std::vector<int> actgroup_idx;

  std::vector<double> grad(d, 0);
  std::vector<double> grad_group(kG, 0);
  std::vector<double> grad_group_master(kG, 0);

  double t=1;

  double a = 0, g = 0, L = 0, sum_r2 = 0;
  double tmp_change = 0, local_change = 0;
  double gamma=0.8;
  r = Y - Xb;


  sum_r2 = r.matrix().dot(r.matrix());
  L = sqrt(sum_r2 / n);

  double dev_thr = fabs(L) * prec;



  for(int kg=0;kg<kG;kg++){



    for(int kk=0; kk<groupLen[kg]; kk++){
      double rXX=(X.col(rangeGroupInd[kg]+kk)*r).sum()/(n*L);
      grad_group[kg]=pow(thresholdl1(rXX, lambda1),2)+grad_group[kg];
    }
    grad_group_master[kg] = grad_group[kg];

  }





  w1_master = w1;
  Xb_master = Xb;


  std::vector<double> stage_group_lambdas(kG, 0);
  std::vector<double> betaIsZero(kG,0);


  double thresh=0.01;




  for(int i=0;i<nlambda;i++)
  {
    w1 = w1_master;
    Xb = Xb_master;


    for(int kg=0;kg<kG;kg++){
      grad_group[kg] = grad_group_master[kg];
      actset_group[kg] = actset_group_master[kg];

    }

    // init the active set
    double threshold;
    if (i > 0)
      threshold = (2 * lambda2[i] - lambda2[i - 1]);
    else
      threshold = (2 * lambda2[i]);

    for (int kg = 0; kg < kG; kg++)
    {
      stage_group_lambdas[kg] = lambda2[i]*sqrt(groupLen[kg]);


      if (grad_group[kg] > threshold) actset_group[kg] = 1;


    }

    r = Y - Xb;

    sum_r2 = r.matrix().dot(r.matrix());
    L = sqrt(sum_r2 / n);


    // loop level 0: multistage convex relaxation
    int loopcnt_level_0 = 0;


    int idx;
    double updated_coord;


    while (loopcnt_level_0 < num_relaxation_round)
    {
      loopcnt_level_0++;

      // loop level 1: active set update
      int loopcnt_level_1 = 0;
      bool terminate_loop_level_1 = true;


      while (loopcnt_level_1 < max_iter)
      {
        loopcnt_level_1++;
        terminate_loop_level_1 = true;


        std::vector<int> isactive_ind(d, 0);
        // initialize actset_idx
        actset_idx.clear();

        actgroup_idx.clear();
        //
        //
        for (int kg = 0; kg < kG; kg++)
        {
          if (actset_group[kg])
          {


            sum_r2 = r.matrix().dot(r.matrix());
            L = sqrt(sum_r2 / n);



            std::vector<double> rXX(groupLen[kg], 0);
            std::vector<double> lXb(groupLen[kg], 0);
            Eigen::VectorXd betatemp;



            betatemp.resize(groupLen[kg]);
            betatemp.setZero();

            std::vector<double> ldot(groupLen[kg], 0);

            for(int kk=0; kk<groupLen[kg]; kk++){
              rXX[kk]=(r * X.col(rangeGroupInd[kg]+kk)).sum()/sqrt(sum_r2);
              ldot[kk]=-(r * X.col(rangeGroupInd[kg]+kk)).sum()/(n*L);

              betatemp[kk]=w1[rangeGroupInd[kg]+kk];

            }


            Eigen::MatrixXd lddot;

            Eigen::MatrixXd xx;
            lddot.resize(groupLen[kg], groupLen[kg]);

            lddot.setZero();

            xx.resize(groupLen[kg], groupLen[kg]);
            xx.setZero();



            for(int kk=0; kk<groupLen[kg]; kk++){
              for(int jj=0; jj<groupLen[kg];jj++){
                xx(kk,jj)=(X.col(rangeGroupInd[kg]+kk)*X.col(rangeGroupInd[kg]+jj)).sum();


                lddot(kk,jj)=(xx(kk,jj)-rXX[kk]*rXX[jj])/(n*L);
              }

              lXb[kk]=(lddot.row(kk)*betatemp).sum();
              grad_group[kg]=grad_group[kg]+pow(thresholdl1(lXb[kk]-ldot[kk],lambda1),2);

            }


            if(grad_group[kg] <= pow(stage_group_lambdas[kg],2))
            {

              betaIsZero[kg] = 1;
              for(int kk=0;kk<groupLen[kg];kk++){
                Xb = Xb - w1[kk + rangeGroupInd[kg]] * X.col(kk+rangeGroupInd[kg]);
              }
              r=Y-Xb;
              for(int kk=0;kk<groupLen[kg];kk++)
              {
                w1[kk + rangeGroupInd[kg]] = 0;
              }
              sum_r2 = 0.0;

              sum_r2 = r.matrix().dot(r.matrix());
              L = sqrt(sum_r2 / n);
            }
            else
            {


              actgroup_idx.push_back(kg);

              Eigen::ArrayXd theta= w1;


              betaIsZero[kg] = 0;
              std::vector<double> z(groupLen[kg], 0);
              std::vector<double> U(groupLen[kg], 0);
              std::vector<double> G(groupLen[kg], 0);
              std::vector<double> lXbeta0(groupLen[kg], 0);


              int count = 0;
              int check = 100000;


              // innergroup
              Eigen::VectorXd beta_t_1;

              beta_t_1.resize(groupLen[kg]);
              beta_t_1.setZero();

              for(int kk=0;kk<groupLen[kg];kk++)
              {
                beta_t_1[kk] = w1[kk + rangeGroupInd[kg]];
              }

              for(int kk=0; kk<groupLen[kg]; kk++){
                lXbeta0[kk]=(lddot.row(kk)*beta_t_1).sum();

              }

              while(count <= innerIter && check > thresh)
              {

                count=count+1;


                Eigen::VectorXd betatemp;
                //

                betatemp.resize(groupLen[kg]);
                betatemp.setZero();

                for(int kk=0;kk<groupLen[kg];kk++){
                  betatemp[kk]=w1[rangeGroupInd[kg]+kk]- beta_t_1[kk];

                }




                std::vector<double> grad(groupLen[kg], 0);

                for(int kk=0; kk<groupLen[kg]; kk++){
                  lXb[kk]=(lddot.row(kk)*betatemp).sum();

                  grad[kk]=lXb[kk]+ldot[kk];

                }


                double diff = -1;


                double l1=0;
                double l2=0;
                for(int kk=0;kk<groupLen[kg];kk++)
                {
                  l1=l1+(ldot[kk]-lXbeta0[kk])*w1[rangeGroupInd[kg]+kk];
                  for(int jj=0;jj<groupLen[kg];jj++){
                    l2=l2+ w1[rangeGroupInd[kg]+kk]*lddot(kk,jj)*w1[rangeGroupInd[kg]+jj];
                  }

                }
                l2=l2/2;

                double Lold =(l1+l2);

                double t=1;
                // each variable inner-group update
                while(diff < 0)
                {
                  for(int kk=0;kk<groupLen[kg];kk++)
                  {


                    z[kk] = w1[kk + rangeGroupInd[kg]] - t * grad[kk];
                    if((z[kk] < lambda1 * t)& (z[kk] > -lambda1* t))
                    {
                      z[kk] = 0;
                    }
                    if(z[kk] > lambda1* t)
                    {
                      z[kk] = z[kk] - lambda1 * t;
                    }
                    if(z[kk] < -lambda1 * t)
                    {
                      z[kk] = z[kk] + lambda1 * t;
                    }
                  }

                  double norm = 0;
                  for(int kk=0;kk<groupLen[kg];kk++)
                  {
                    norm = norm + pow(z[kk],2);
                  }


                  norm = sqrt(norm);
                  double uOp;
                  if(norm != 0){
                    uOp = (1 - stage_group_lambdas[kg]*t/norm);
                  }
                  else{uOp = 0;}

                  if(uOp < 0)
                  {
                    uOp = 0;
                  }

                  for(int kk=0;kk<groupLen[kg];kk++)
                  {
                    U[kk] = uOp*z[kk];
                    G[kk] = 1/t *(w1[kk + rangeGroupInd[kg]] - U[kk]);


                  }

                  double l1=0;
                  double l2=0;
                  for(int kk=0;kk<groupLen[kg];kk++)
                  {
                    l1=l1+(ldot[kk]-lXbeta0[kk])*U[kk];
                    for(int jj=0;jj<groupLen[kg];jj++){
                      l2=l2+ U[kk]*lddot(kk,jj)*U[jj];
                    }
                  }
                  l2=l2/2;

                  double Lnew =(l1+l2);

                  double sqNormG = 0;
                  double iProd = 0;

                  for(int kk=0;kk<groupLen[kg];kk++)
                  {
                    sqNormG = sqNormG + pow(G[kk],2);
                    iProd = iProd + grad[kk] * G[kk];
                  }

                  diff = Lold - Lnew - t * iProd + t/2 * sqNormG;

                  t = t * gamma;

                }

                double check = 0;
                int reset=10;
                for(int kk=0;kk<groupLen[kg];kk++)
                {

                  check = check + fabs(theta[kk + rangeGroupInd[kg]] - U[kk]);
                  double hh= double(count%reset)/double(count%reset+3);
                  w1[kk + rangeGroupInd[kg]] = theta[kk + rangeGroupInd[kg]] + hh * (U[kk] - theta[kk + rangeGroupInd[kg]]);
                  theta[kk + rangeGroupInd[kg]] = U[kk];


                }

              }


              std::vector<double> tmp(groupLen[kg], 0);

              for(int kk=0;kk<groupLen[kg];kk++)
              {
                tmp[kk]= w1[kk + rangeGroupInd[kg]]-beta_t_1[kk];
              }
              for(int kk=0;kk<groupLen[kg];kk++){
                Xb = Xb + tmp[kk] * X.col(kk+rangeGroupInd[kg]);
                r=r-tmp[kk]*X.col(kk+rangeGroupInd[kg]);
                updated_coord = w1[kk + rangeGroupInd[kg]];
                if (fabs(updated_coord) > 0) {
                  actset_idx.push_back(kk + rangeGroupInd[kg]);
                  isactive_ind[kk + rangeGroupInd[kg]]=1;
                }
              }

              sum_r2 = 0.0;

              sum_r2 = r.matrix().dot(r.matrix());
              L = sqrt(sum_r2 / n);
            }

          }
        }



        // loop level 2: proximal newton on active set
        int loopcnt_level_2 = 0;

        bool terminate_loop_level_2 = true;

        while (loopcnt_level_2 < max_iter)
        {
          loopcnt_level_2++;
          terminate_loop_level_2 = true;

          for (unsigned int k = 0; k < actgroup_idx.size(); k++)
          {

            //update the ldot and lddot-----due to 21-11-9
            double kg = actgroup_idx[k];

            //cout<<kg<<endl;
            sum_r2 = r.matrix().dot(r.matrix());
            L = sqrt(sum_r2 / n);


            Eigen::ArrayXd theta= w1;

            std::vector<double> z(groupLen[kg], 0);
            std::vector<double> U(groupLen[kg], 0);
            std::vector<double> G(groupLen[kg], 0);
            std::vector<double> lXbeta0(groupLen[kg], 0);



            std::vector<double> rXX(groupLen[kg], 0);


            std::vector<double> ldot(groupLen[kg], 0);

            for(int kk=0; kk<groupLen[kg]; kk++){
              rXX[kk]=(r * X.col(rangeGroupInd[kg]+kk)).sum()/sqrt(sum_r2);
              ldot[kk]=-(r * X.col(rangeGroupInd[kg]+kk)).sum()/(n*L);

            }
            Eigen::MatrixXd lddot;
            Eigen::MatrixXd xx;
            lddot.resize(groupLen[kg], groupLen[kg]);

            lddot.setZero();

            xx.resize(groupLen[kg], groupLen[kg]);
            xx.setZero();


            for(int kk=0; kk<groupLen[kg]; kk++){
              for(int jj=0; jj<groupLen[kg]; jj++){
                xx(kk,jj)=(X.col(rangeGroupInd[kg]+kk)*X.col(rangeGroupInd[kg]+jj)).sum();
                lddot(kk,jj)=(xx(kk,jj)-rXX[kk]*rXX[jj])/(n*L);
              }

            }


            int count = 0;
            int check = 100000;


            // innergroup
            Eigen::VectorXd beta_t_1;

            beta_t_1.resize(groupLen[kg]);
            beta_t_1.setZero();

            for(int kk=0;kk<groupLen[kg];kk++)
            {

              beta_t_1[kk] = w1[kk + rangeGroupInd[kg]];

            }

            while(count <= innerIter && check > thresh)
            {
              count=count+1;

              std::vector<double> grad(groupLen[kg], 0);
              std::vector<double> lXbeta0(groupLen[kg], 0);
              std::vector<double> lXb(groupLen[kg], 0);


              Eigen::VectorXd betatemp;
              betatemp.resize(groupLen[kg]);
              betatemp.setZero();

              for(int kk=0; kk<groupLen[kg]; kk++){
                if(isactive_ind[kk+ rangeGroupInd[kg]]){
                  betatemp[kk]=w1[rangeGroupInd[kg]+kk]- beta_t_1[kk];

                  lXbeta0[kk]=(lddot.row(kk)*beta_t_1).sum();

                }
              }



              for(int kk=0; kk<groupLen[kg]; kk++){
                if(isactive_ind[kk+ rangeGroupInd[kg]]){
                  lXb[kk]=(lddot.row(kk)*betatemp).sum();
                  grad[kk]=lXb[kk]+ldot[kk];
                }
              }



              double diff = -1;



              double l1=0;
              double l2=0;
              for(int kk=0;kk<groupLen[kg];kk++)
              {
                l1=l1+(ldot[kk]-lXbeta0[kk])*w1[rangeGroupInd[kg]+kk];
                for(int jj=0;jj<groupLen[kg];jj++){
                  l2=l2+ w1[rangeGroupInd[kg]+kk]*lddot(kk,jj)*w1[rangeGroupInd[kg]+jj];
                }

              }
              l2=l2/2;

              double Lold =(l1+l2);


              std::vector<double> z(groupLen[kg], 0);
              std::vector<double> U(groupLen[kg], 0);
              std::vector<double> G(groupLen[kg], 0);




              // each variable inner-group update
              double t=1;
              while(diff < 0)
              {
                for(int kk=0;kk<groupLen[kg];kk++)
                {

                  if(isactive_ind[kk+ rangeGroupInd[kg]]){
                    z[kk] = w1[kk+ rangeGroupInd[kg]] - t * grad[kk];
                    if((z[kk] < lambda1 * t) &(z[kk] > -lambda1* t))
                    {
                      z[kk] = 0;
                    }
                    if(z[kk] > lambda1* t)
                    {
                      z[kk] = z[kk] - lambda1 * t;
                    }
                    if(z[kk] < -lambda1 * t)
                    {
                      z[kk] = z[kk] + lambda1 * t;
                    }
                  }else{
                    z[kk]=0;
                  }
                }

                double norm = 0;
                for(int kk=0;kk<groupLen[kg];kk++)
                {
                  norm = norm + pow(z[kk],2);
                }
                norm = sqrt(norm);
                double uOp=0;
                if(norm != 0){
                  uOp = (1 - stage_group_lambdas[kg]*t/norm);
                }
                else{uOp = 0;}

                if(uOp < 0)
                {
                  uOp = 0;
                }

                for(int kk=0;kk<groupLen[kg];kk++)
                {
                  U[kk] = uOp*z[kk];
                  G[kk] = 1/t *(w1[kk+ rangeGroupInd[kg]] - U[kk]);

                }

                double l1=0;
                double l2=0;
                for(int kk=0;kk<groupLen[kg];kk++)
                {
                  l1=l1+(ldot[kk]-lXbeta0[kk])*U[kk];
                  for(int jj=0;jj<groupLen[kg];jj++){
                    l2=l2+ U[kk]*lddot(kk,jj)*U[jj];
                  }

                }
                l2=l2/2;

                double Lnew =(l1+l2);

                double sqNormG = 0;
                double iProd = 0;

                for(int kk=0;kk<groupLen[kg];kk++)
                {
                  sqNormG = sqNormG + pow(G[kk],2);
                  iProd = iProd + grad[kk] * G[kk];
                }

                diff = Lold - Lnew - t * iProd + t/2 * sqNormG;

                t = t * gamma;
              }



              double check = 0;

              int reset=10;
              for(int kk=0;kk<groupLen[kg];kk++)
              {

                if(isactive_ind[kk+ rangeGroupInd[kg]]){
                  check = check + fabs(theta[kk + rangeGroupInd[kg]] - U[kk]);
                  double hh= double(count%reset)/double(count%reset+3);
                  w1[kk + rangeGroupInd[kg]] =  theta[kk + rangeGroupInd[kg]] + hh*(U[kk] - theta[kk + rangeGroupInd[kg]]);
                  theta[kk + rangeGroupInd[kg]] = U[kk];

                }
              }
            }


            std::vector<double> tmp(groupLen[kg], 0);
            for(int kk=0;kk<groupLen[kg];kk++)
            {
              tmp[kk]= w1[kk + rangeGroupInd[kg]]-beta_t_1[kk];
            }


            for(int kk=0;kk<groupLen[kg];kk++){
              Xb = Xb + tmp[kk] * X.col(kk+rangeGroupInd[kg]);
              r=r-tmp[kk]*X.col(kk+rangeGroupInd[kg]);


            }

            sum_r2 = 0.0;

            sum_r2 = r.matrix().dot(r.matrix());
            L = sqrt(sum_r2 / n);



            double tmp_change=0;
            double local_change=0;
            for(int kk=0;kk<groupLen[kg];kk++){
              if(isactive_ind[kk+ rangeGroupInd[kg]]){

                tmp_change = beta_t_1[kk]- w1[kk+ rangeGroupInd[kg]];
                double a =  (X.col(kk+ rangeGroupInd[kg]) * X.col(kk+ rangeGroupInd[kg]) * (1 - r * r/(L*L*n))).sum()/(n*L);
                local_change = a * tmp_change * tmp_change / (2 * L * n)+local_change;

              }
            }
            if (local_change > dev_thr){
              terminate_loop_level_2 = false;
            }
          }


          if (terminate_loop_level_2)
            break;

        }
        bool new_active_idx=false;
        // check stopping criterion 2: active set change
        bool new_group_idx = false;
        for (int kg = 0; kg < kG; kg++)
          if (actset_group[kg] == 0)
          {
            grad_group[kg]=0;
            std::vector<double> rXX(groupLen[kg], 0) ;

            for(int kk=0; kk<groupLen[kg]; kk++){
              rXX[kk]=(X.col(rangeGroupInd[kk]+kk)*r).sum()/(n*L);
              grad_group[kg]=pow((thresholdl1(rXX[kk], lambda1)),2)+grad_group[kg];
            }

            gr[kg] = fabs(grad_group[kg]);
            bool new_active_idx;
            if (gr[kg] > stage_group_lambdas[kg])
            {
              actset_group[kg] = 1;
              new_active_idx = true;
            }
          }

          if(!new_active_idx)
            break;



      }

      if (loopcnt_level_0 == 1)
      {
        for (int j = 0; j < d; j++)
        {
          w1_master[j] = w1[j];



        }
        for(int kg=0;kg<kG;kg++){
          grad_group_master[kg] = grad_group[kg];
          actset_group_master[kg] = actset_group[kg];
        }
        for (int j = 0; j < n; j++) Xb_master[j] = Xb[j];
      }

    }

  }


  if(standard==1){
    for(int ii=0;ii<d;ii++){

      w1[ii]=w1[ii]/varX[ii];
    }
  }


  return w1;
}


double thresholdl1(double x, double thr) {
  if (x > thr)
    return x - thr;
  else if (x < -thr)
    return x + thr;
  else
    return 0;
}
//data: original analized matrix of nrow X p
//K: the group size; tmp: weight matrix; index: group vector;
//lambda1: lasso tuning; lambda2: group tuning;
//nlam: sqrt_regression step's group parameter's number, also the length of lambda2;
//groupLen: group length; rangeGroupInd:
//beta: the initialized value of ordered coef;
//f: ordered function
//tau: cutoff
//cc: the intercepts K*p
//[[Rcpp::export]]
List nodewise(Eigen::Map<Eigen::MatrixXd> data, int K, int nrow, int p,int num_relaxation_round,
              Eigen::Map<Eigen::MatrixXd> residual,
              int max_iter,
              int innerIter,
              Eigen::Map<Eigen::MatrixXd> tmp,
              NumericVector index,
              NumericVector lambda1,
              NumericVector lambda2,
              int nlam,
              NumericVector groupLen,
              NumericVector rangeGroupInd,
              Eigen::Map<Eigen::MatrixXd> beta,
              Function f,
              double tau,
              int standard,
              Eigen::Map<Eigen::MatrixXd> cc)
{

  Eigen::ArrayXd out;
  out.resize(K*(p-1));

  Eigen::ArrayXd outcut;
  outcut.resize(K*(p-1));


  Eigen::ArrayXXd betatmp=beta;

  Eigen::MatrixXd  weightX;

  weightX.resize(nrow, (K*p));
  weightX.setZero();


  List outt(p);
  List outt2(p);

  Eigen::MatrixXd  XX; // the predictor matrix of dimension n*KX(p-1)
  XX.resize((nrow*K),K*(p-1));
  XX.setZero();

  Eigen::MatrixXd  Y;// the predicted matrix of dimension n*KX(1)
  Y.resize((nrow*K),1);
  Y.setZero();



  int l=0;
  for(int k=0;k<K;k++){
    for(int j=0;j<p;j++){
      weightX.col(l)=tmp.col(k).cwiseProduct(data.col(j));

      l++;
    }
  }
  for(int i=0;i<(p);i++){


    for(int k=0;k<K;k++){
      if(i<1){

        XX.block((k*nrow), (k*(p-1)), (nrow), (p-1))=weightX.block(0, (k*p+1), (nrow), (p-1))/residual(k,i);

        Y.block((k*nrow), 0, (nrow), 1)=(weightX.col(k*p)-cc(k,i)*tmp.col(k))/residual(k,i);
      }else{
        XX.block((k*nrow), (k*(p-1)), (nrow), (i))=weightX.block(0, (k*(p)), (nrow), (i))/residual(k,i);
        XX.block((k*nrow), (k*(p-1)+i), (nrow), (p-i-1))=weightX.block(0, (k*p+i+1), (nrow), (p-i-1))/residual(k,i);
        Y.block((k*nrow), 0, (nrow), 1)=(weightX.col(k*p+i)-cc(k,i)*tmp.col(k))/residual(k,i);
      }
    }

    Eigen::ArrayXd btmp;
    btmp=betatmp.col(i);

    NumericMatrix XX1;


    XX1=f(index,XX);

    Eigen::MatrixXd XX2;
    XX2=Test(XX1);



    Eigen::MatrixXd onesM;
    onesM.resize((nrow*K),1);






    out=sqrt_regression(btmp,XX2,Y,groupLen,rangeGroupInd,lambda2,lambda1[0],nlam,(K*(p-1)),(p-1),num_relaxation_round,max_iter,innerIter,standard);

    outcut=out;

    for(int k = 0; k < (K*(p-1)); k++)
      if(fabs(outcut[k])<tau){
        outcut[k]= 0;
      }

      Eigen::MatrixXd temp;
      temp.resize((nrow*K), 1);
      for(int j = 0; j < (nrow*K); j++)
      {
        temp(j, 0) = 0;
        for(int k = 0; k < (K*(p-1)); k++)
          temp(j, 0) += XX2(j, k)*outcut[k];
        temp(j, 0) = Y(j,0) - temp(j, 0);
      }


      outt[i]=out;
      outt2[i]=temp;

  }
  return List::create(
    _["coeflist"]=outt,
    _["residual"]=outt2
  );
}



Eigen::MatrixXd Test(NumericMatrix AA){

  Map<Eigen::MatrixXd> A(as<Map<Eigen::MatrixXd> >(AA));

  return(A);
}
