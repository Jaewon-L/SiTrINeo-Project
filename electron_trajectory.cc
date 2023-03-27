
const int L = 20; // magnetic field length(mm)
const double B = 0.5; // magnatic field strength(T)
const double angle = 0.5; // angle range setting(rad)

double Array_max(double a[]) {
  double max = a[0];
  for(int b = 1; b <= L; b++) {
    if(a[b] > max) {
      max = a[b];
    }
  }
  return max;
}

void electron_trajectory() {

  const int SIZE = 275; // momentum size
  const int AS = 121; // angle size
  const double Mev = 0.511; // mass of electron(MeV/c^2)
  const double Mek = 9.109 * pow(10, -31); // mass of electron(kg)
  const double Q = 1.602 * pow(10, -19); // charge of electron(C)
  double V[SIZE]; // electron of velocity(c)
  double VM[SIZE];// electron of velocity(m/s)
  double R[SIZE]; // Radius(m)
  double P[SIZE]; // momentum of electron(MeV/c)
  double alpha[SIZE]; // angle of electron entering sensor(radian)
  double beta = 0;

  int i = 0;

  for(double j = 0; j < 2.74; j += 0.01) {
    P[i] = j;
    i++;
  }

  for(int k = 0; k < SIZE; k++){
    V[k] = P[k] / sqrt(pow(P[k], 2) + pow(Mev, 2));
   // printf("%d: %f: %f \n", k, P[k], V[k]);
    VM[k] = V[k] * 2.998 * pow(10, 8);
  }

  TCanvas * c1 = new TCanvas("c1", "multipads", 1500, 1000);

  c1->Divide(2,2);

  c1->cd(1);

  TGraph * grp = new TGraph(); 
  grp->SetTitle("Momentum-Radius relation; Momentum(MeV); Radius(mm)"); 

  for(int l = 0; l < SIZE; l++) {
    R[l] = (Mek * VM[l]) / (Q * B);
   // printf("%d: P: %fMeV, V: %fc, R:%fmm \n", l, P[l], V[l], R[l] * 1000);
    grp->SetPoint(l, P[l], R[l] * 1000);  
  }
  grp->Draw("APL");

// angle range setting
  const double IRP = atan(angle) / 120; // 27/60'
  //const double IRM = atan(-0.3) / 120; // 27/60'
  alpha[0] = 0;

  for(int m = 1; m < AS; m++) {
    alpha[m] = alpha[m-1] + IRP;
  }

  int s = 0;
  double FY; 
  double SY[L+1];
  double Y;
  double S_2Y;
  double S_1Y;
  double Y1;
  double Y2;

//initial version
//----------------------------------------------------------------------

   c1->cd(2);  

   TH1D * value1 = new TH1D("value1","Y_coord theoretical sensor1_first",150,-150,30);


   for(int n = 0; n < SIZE; n++) {
     for(int o = 0; o < AS; o++) {
       beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
       Y1 = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2) + 10 * tan(beta);
       //printf("n: %d, o: %d, s: %d: Y: %fmm, P: %fMeV, alpha: %frad \n",n ,o, s, Y, P[n], alpha[o]);
      // s++;
       value1->SetXTitle("Y_coord_sensor1(mm)");
       value1->SetYTitle("event");
       value1->Fill(Y1);
     }
   }
   value1->Draw();
   c1->Draw(); 

//-------------------------------------------------------------------
//second version

   c1->cd(3);  

   TH1D * value2 = new TH1D("value2","Y_coord theoretical sensor1_second",150,-150,30);

   for(int n = 0; n < SIZE; n++) {
     for(int o = 0; o < AS; o++) {
       beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
       FY = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2);
       if(isnan(FY) != 0) {
         continue;
       }
       else { 
	 for(int q = 0; q <= L; q++) {      
           SY[q] = 20 * tan(alpha[o]) + q * tan((alpha[o] + beta) / 2);
         }  
         if(Array_max(SY) > 10) {
	   continue;
	 }
         else {
           beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
           Y2 = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2) + 5 * tan(beta);
           value2->Fill(Y2);
         }
       }
       //printf("n: %d, o: %d, s: %d: Y: %fmm, P: %fMeV, alpha: %frad \n",n ,o, s, Y[s], P[n], alpha[o]);    
     }
   }      
   value2->SetXTitle("Y_coord_sensor1(mm)");
   value2->SetYTitle("event");
   value2->Draw();
   c1->Draw(); 

//------------------------------------------------------------------------

//final version

   c1->cd(4);  

   TH1D * value3 = new TH1D("value3","Y_coord theoretical sensor1_final",150,-150,30);
   FILE * fp = fopen("senor_value.txt", "w"); 
  
   for(int n = 0; n < SIZE; n++) {
     for(int o = 0; o < AS; o++) {
       beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
       FY = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2);
       if(isnan(FY) != 0) {
         continue;
       }
       else { 
	 for(int q = 0; q <= L; q++) {      
           SY[q] = 20 * tan(alpha[o]) + q * tan((alpha[o] + beta) / 2);
         }  
         if(Array_max(SY) > 10) {
	   continue;
	 }
         else {
           beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
           S_2Y = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2) + 5 * tan(beta);
           if(-20 < S_2Y && S_2Y < 0) {
             beta = asin(sin(alpha[o]) - ((0.3 * B * L) / P[n]));
             S_1Y = 20 * tan(alpha[o]) + L * tan((alpha[o] + beta) / 2) + 10 * tan(beta);
             if(-20 < S_1Y && S_1Y < 0) {
               printf("n: %d, o: %d, s: %d: Y: %fmm, P: %fMeV, alpha: %frad \n",n ,o, s, S_1Y, P[n], alpha[o]);
               value3->Fill(S_1Y);
               fprintf(fp, "n: %d, o: %d, s: %d: Y: %fmm, P: %fMeV, alpha: %frad \n",n ,o, s, S_1Y, P[n], alpha[o]);
             }
           } 
         }
       }
     s++;  //printf("n: %d, o: %d, s: %d: Y: %fmm, P: %fMeV, alpha: %frad \n",n ,o, s, Y[s], P[n], alpha[o]);    
     }
   }      
   fclose(fp);
   value3->SetXTitle("Y_coord_sensor1(mm)");
   value3->SetYTitle("event");
   value3->Draw();
   c1->cd(2);
   value3->Draw("same");

  // c1->Draw(); 
}












