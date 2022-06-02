double *Beasley_Springer_Moro(int n, double u[]){
    double y[n];
    double r[n];
    double *x = new double[n];
    // Constants needed for algo
    double a_0 = 2.50662823884;     double b_0 = -8.47351093090;
    double a_1 = -18.61500062529;   double b_1 = 23.08336743743;
    double a_2 = 41.39119773534;    double b_2 = -21.06224101826;
    double a_3 = -25.44106049637;   double b_3 = 3.13082909833;

    double c_0 = 0.3374754822726147; double c_5 = 0.0003951896511919;
    double c_1 = 0.9761690190917186; double c_6 = 0.0000321767881768;
    double c_2 = 0.1607979714918209; double c_7 = 0.0000002888167364;
    double c_3 = 0.0276438810333863; double c_8 = 0.0000003960315187;
    double c_4 = 0.0038405729373609;

    for(int i = 0; i < n; i++){
        y[i] = u[i] - 0.5;
    }
    for(int i = 0; i < n; i++){
        if(fabs(y[i]) < 0.42){
            r[i] = y[i]*y[i];
            x[i] = y[i]*(((a_3 * r[i] + a_2)*r[i] + a_1)*r[i] + a_0)/((((b_3 * r[i] + b_2)*r[i] + b_1)*r[i] + b_0)*r[i] + 1.0);
        } else{
            r[i] = u[i];
            if(y[i] > 0){
                r[i] = 1 - u[i];
            }
            r[i] = log(-log(r[i]));
            x[i] = c_0 + r[i]*(c_1 + r[i]*(c_2 + r[i] *(c_3 + r[i] * (c_4 + r[i] * (c_5 + r[i]*(c_6 + r[i] * (c_7 + r[i]*c_8)))))));
            if(y[i] < 0){
                x[i] = -x[i];
            }
        }
    }
    ofstream myfile ("bb.txt");
    for(int j = 0; j < n; j++){
        myfile << x[j] << ";" << endl;
    }
    return x;
}
