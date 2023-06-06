
    printf("\n\n\n*****************************************************\n");
    printf("*    Changing directory and creating files for      *\n");
    printf("*                     the output                    *\n");
    printf("*****************************************************\n");
    ofstream fp_Loss, fp_E;
    multivec1D<ofstream> fp_TAk; 
    ofstream fp_J1, fp_J2, fp_J3, fp_TAbs, fp_Jintra; 
    ofstream fp_wf;
    

    if (rank_==0)
    {
        // system("mkdir Output");
        fp_Loss.open("Output/Losses.txt");
        printf("->     Losses.txt opened  \n");
         
        fp_E.open("Output/EF.txt");
        printf("->     EF.txt opened  \n");


        Coulomb_set.exciton_file.open("Output/Exciton.txt");
        
        if(iCurrent)
        {
            fp_J1.open("Output/J1.txt");
            printf("->     J1.txt opened  \n");
    
            fp_J2.open("Output/J2.txt");
            printf("->     J2.txt opened  \n");
        
            fp_J3.open("Output/J3.txt");
            printf("->     J3.txt opened  \n");
            
            fp_Jintra.open("Output/Jintra.txt");
            printf("->     Jintra.txt opened  \n");
        }
        if(iTAbs)
        {
            fp_TAbs.open("Output/TransientAbs.txt");
            printf("->     TransientAbs.txt opened  \n");
        }
        if(iTAbsK)
        {
            if(nTAk!=0) fp_TAk.resize(nTAk);
            for(int i=0; i<nTAk; i++)
            {
                stringstream sname;
                sname.seekp(0,ios::beg); sname << TAkpt[i][0] << "_" << TAkpt[i][1] << "_" << TAkpt[i][2];       
                string name_file1= "Output/TAbsK_" + sname.str() + ".txt"; 
                string name_file2= "TAbsK_" + sname.str() + ".txt"; 
               fp_TAk[i].open(name_file1.c_str());
                printf("->     %s opened  \n", name_file2.c_str());
    
            }
        }

        if (Coulomb_set.Print_new_band_dispersion or Diff_Eq.PrintPopulation){
            
            string name_file; 
            name_file= "Output/t_Ek.txt";
            t_Ek.open(name_file.c_str());
        }

    }