// Metod to solve a set of differentional equations with standart Runge - Kutta 4

for(int it=iti;it<=int(itfi/4); it++)
  //for(int it=iti;it<=(iti+20);it++)
  {
      time=it*dt;
      // if resolution condition is true, we save some parameters here
      if (it % (1*it_resolution) == 0)
      {
        for (int ik = begin_count; ik < end_count; ik++){ // all local wave vectors
          for (int ic=0; ic<Ncv; ic++){ //  over bands
            #pragma omp simd
            for (int jc=0; jc<Ncv; jc++){ //
              P0[ik][ic][jc] = P0_local[ik - begin_count][ic][jc];
            }
          }
        }
        #pragma omp barrier // sinchronise threads
         #include "if_resolution.cpp"
      }  // if resolution condition closed


      if (id_masters[thr_total - 4] == true){
       // this is shared variables, only one thread calculate it
        EF1 = pulse1.E(time);
        EF2 = pulse2.E(time);
      }



      //------ 1st step Runge-Kutta
      #pragma omp barrier // sinchronise threads
      if(iTightBinding) Runge_Kutta_Df(TB_Model,k_local,P0_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      else              Runge_Kutta_Df(WModel,  k_local,P0_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ad(P0_local,P1_local,Pv_local,dt6);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ad(P0_local,P2_local,Pv_local,dt2);


      //------ 2nd step Runge-Kutta
      time=it*dt + dt2;
      if (id_masters[thr_total - 4] == true)
      { // this is shared variables
      EF1 = pulse1.E(time);
      AF1 = pulse1.A_crys(time);
      EF2 = pulse2.E(time);
      }

      #pragma omp barrier // sinchronise threads
      if(iTightBinding) Runge_Kutta_Df(TB_Model,k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      else              Runge_Kutta_Df(WModel,k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ac(P1_local,Pv_local,dt3);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ad(P0_local,P2_local,Pv_local,dt2);


      #pragma omp barrier // sinchronise threads
      //------ 3rd step Runge-Kutta
      if(iTightBinding) Runge_Kutta_Df(TB_Model,k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      else              Runge_Kutta_Df(WModel,  k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ac(P1_local, Pv_local, dt3);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ad(P0_local, P2_local, Pv_local, dt);

      //------ 4th step Runge-Kutta
      time=(it+1)*dt;
      if (id_masters[thr_total - 4] == true)
      { // this is shared variables
      EF1 = pulse1.E(time);
      AF1 = pulse1.A_crys(time);
      EF2 = pulse2.E(time);
      }

      #pragma omp barrier // sinchronise threads
      if(iTightBinding) Runge_Kutta_Df(TB_Model,k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      else              Runge_Kutta_Df(WModel,  k_local,P2_local,Pv_local,T,Nb,EF1,EF2,pulse2.wl,AF1,Coulomb_set, id_masters, integrWeight_local);
      #pragma omp barrier // sinchronise threads
      Runge_Kutta_Ad(P1_local, P0_local, Pv_local,dt6);
      #pragma omp barrier // sinchronise threads






    }//end time evolution