// time loop with constant dt

for(int it_pr=iti;it_pr<=itfi;it_pr++)
//for(int it=iti;it<= 100;it++)
// int it = iti;
{   
    #pragma omp master
    { // shared variable, only master change it
        time_loop=it_pr*dt;
        it = it_pr;
    }
    #pragma omp barrier

    if (it % it_resolution == 0){   
        #include "if_resolution.cpp" // printig all we need
    } // end if it_resolution
    
    if(Diff_Eq.Taylor){
        #include "Taylor_DE_solver.cpp"
    } else{
        // RK time evolution
        #include "RungeKutta.cpp" 
    }
    

}//end time evolution