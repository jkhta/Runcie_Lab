#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>


// [[Rcpp::export]]
List c_jaeger_model_V3(NumericVector t, NumericVector X, List parms_list) {
	//NumericVector X,
	// parms_list) 
	//double K_13           = parms_list["K_13"];
        NumericVector times_model   = t;
        NumericVector protein_conc  = X;
        double K_13     = parms_list["K_13"];
        double K_23     = parms_list["K_23"];
        double K_4_3    = parms_list["K_4_3"];
        double K_23_4   = parms_list["K_23_4"];
        double K_13_4   = parms_list["K_13_4"];
        double K_23_5   = parms_list["K_23_5"];
        double K_13_5   = parms_list["K_13_5"];
        double K_4_5    = parms_list["K_4_5"];
        double K_5_4    = parms_list["K_5_4"];
        double h_4_3    = parms_list["h_4_3"];
        double h_23_4   = parms_list["h_23_4"];
        double h_13_4   = parms_list["h_13_4"];
        double h_23_5   = parms_list["h_23_5"];
        double h_13_5   = parms_list["h_13_5"];
        double h_4_5    = parms_list["h_4_5"];
        double h_5_4    = parms_list["h_5_4"];
        double h_5_2    = parms_list["h_5_2"];
        double eta_leaf = parms_list["eta_leaf"];
        double eta_leaf2 = parms_list["eta_leaf2"];
        double T_f      = parms_list["T_f"];
        
        NumericVector delta   = as<NumericVector>(parms_list["delta"]); 
        NumericVector v_35S   = as<NumericVector>(parms_list["v_35S"]); 
        NumericVector v1      = as<NumericVector>(parms_list["v1"]); 
        NumericVector v2      = as<NumericVector>(parms_list["v2"]); 
        NumericVector v3      = as<NumericVector>(parms_list["v3"]); 
        NumericVector mutants = as<NumericVector>(parms_list["mutants"]); 
    
    //Hub protein-Protein Binding
    //FT:FD
        double x_13 = K_23*X[0]*X[2]/(K_13*K_23 + K_13*X[1] + K_23*X[0]);
    //TFL1:FD
        double x_23 = K_13*X[1]*X[2]/(K_13*K_23 + K_13*X[1] + K_23*X[0]);
    
	//Hub Gene Activation
    //LFY -> FD
        double p_4_3 = pow(X[3],h_4_3) / (pow(K_4_3,h_4_3) + pow(X[3],h_4_3));
    //FT:FD -> LFY
        double p_13_4 = pow(K_23_4,h_23_4) * pow(x_13,h_13_4) / (pow(K_13_4,h_13_4) * pow(K_23_4,h_23_4) + pow(K_23_4,h_23_4)*pow(x_13,h_13_4) + pow(K_13_4,h_13_4)*pow(x_23,h_23_4));
    //TFL1:FD -> LFY
        double p_23_4 = pow(K_13_4,h_13_4) * pow(x_23,h_23_4) / (pow(K_13_4,h_13_4) * pow(K_23_4,h_23_4) + pow(K_23_4,h_23_4)*pow(x_13,h_13_4) + pow(K_13_4,h_13_4)*pow(x_23,h_23_4));
    //AP1 -> LFY
        double p_5_4 = pow(X[4],h_5_4) / (pow(K_5_4,h_5_4) + pow(X[4],h_5_4));
    //FT:FD -> AP1
        double p_13_5 = pow(K_23_5,h_23_5) * pow(x_13,h_13_5) / (pow(K_13_5,h_13_5) * pow(K_23_5,h_23_5) + pow(K_23_5,h_23_5)*pow(x_13,h_13_5) + pow(K_13_5,h_13_5)*pow(x_23,h_23_5));
    //TFL1:FD -> AP1
        double p_23_5 = pow(K_13_5,h_13_5) * pow(x_23,h_23_5) / (pow(K_13_5,h_13_5) * pow(K_23_5,h_23_5) + pow(K_23_5,h_23_5)*pow(x_13,h_13_5) + pow(K_13_5,h_13_5)*pow(x_23,h_23_5));
    //LFY -> AP1
        double p_4_5 = pow(X[3],h_4_5) / (pow(K_4_5,h_4_5) + pow(X[3],h_4_5));
    
	//Synthesis Rates
	  NumericMatrix rho(5,3);
    
	  rho(0,0) = 1;
	  rho(0,1) = 0;
	  rho(0,2) = 0;
    
	  rho(1,1) = pow(T_f,h_5_2) / (pow(T_f,h_5_2)  + pow(X[4],h_5_2));
	  rho(1,2) = 0;
	  rho(1,0) = 1; // - rho(1,1); Correction not done in Jaeger et al.
    
	  rho(2,1) = p_4_3;
	  rho(2,0) = 1 - rho(2,1);
	  rho(2,2) = 0;
    
	  rho(3,0) = (1 - p_13_4 - p_23_4) * (1-p_5_4);
	  rho(3,1) = p_13_4*(1-p_5_4) + (1 - p_13_4 - p_23_4)*p_5_4;
	  rho(3,2) = p_13_4*p_5_4;
    
	  rho(4,0) = (1 - p_13_5 - p_23_5) * (1-p_4_5);
	  rho(4,1) = p_13_5*(1-p_4_5) + (1 - p_13_5 - p_23_5)*p_4_5;
	  rho(4,2) = p_13_5*p_4_5;
    
	// FT not calculated here
	// eta_leaf = 0.01
	// v1[0] = eta_leaf * t;
    
	  NumericVector Vout(7);
    
	  Vout[0] = v_35S[0] + (rho(0,0)*eta_leaf*X[5] + rho(0,1)*v2[0] + rho(0,2)*v3[0]) * (1-mutants[0]) - (delta[0] * X[0]);
	  Vout[1] = v_35S[1] + (rho(1,0)*v1[1] + rho(1,1)*v2[1] + rho(1,2)*v3[1]) * (1-mutants[1]) - (delta[1] * X[1]);
	  Vout[2] = v_35S[2] + (rho(2,0)*v1[2] + rho(2,1)*v2[2] + rho(2,2)*v3[2]) * (1-mutants[2]) - (delta[2] * X[2]);
	  Vout[3] = v_35S[3] + (rho(3,0)*v1[3] + rho(3,1)*v2[3] + rho(3,2)*v3[3]) * (1-mutants[3]) - (delta[3] * X[3]);
	  Vout[4] = v_35S[4] + (rho(4,0)*v1[4] + rho(4,1)*v2[4] + rho(4,2)*v3[4]) * (1-mutants[4]) - (delta[4] * X[4]);
	  Vout[5] = X[6];
	  
	  // Rcout << v_35S[0] << " ";
	  // Rcout << rho(0,0) << " ";
	  // Rcout << v1[0] << " ";
	  // Rcout << rho(0,1) << " ";
	  // Rcout << v2[0] << " ";
	  // Rcout << rho(0,2) << " ";
	  // Rcout << v3[0] << "\n";
	  // 
	  // 
	  // Rcout << Vout[0] << " ";
	  // Rcout << Vout[1] << " ";
	  // Rcout << Vout[2] << " ";
	  // Rcout << Vout[3] << " ";
	  // Rcout << Vout[4] << " ";
	  // Rcout << Vout[5] << " ";
	  // Rcout << Vout[6] << "\n\n";
	  
        //NumericVector protein_conc  = X;
    return List::create(_["Derivative"]=Vout, _["States"] = NumericVector::create(_["x_13"] = x_13, _["x_23"] = x_23, _["p_13_4"] = p_13_4, _["p_13_5"] = p_13_5, _["p_23_4"] = p_23_4, _["p_23_5"] = p_23_5, _["p_4_3"] = p_4_3, _["p_4_5"] = p_4_5, _["p_5_4"] = p_5_4));
	        //_["parm"] = K_13,
	        //_["states"] = times_model);
	        //_["death"] = protein_conc);
}