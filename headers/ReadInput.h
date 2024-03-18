void Separate_string(string& s, vector<string>& str, size_t& pos)
{
  str.clear();
        s = s.substr(0,pos); //cout << s << endl;
        str.push_back("");
        size_t count = 0;
        size_t w=0;  
        while(w != pos)            
        {
            if( s[w] != ' ')
            {
                str[count] += s[w];
            }
            else
            {
              if( 0 != strcasecmp(str[count].c_str(), "") )
              {
                count++;
                str.push_back("");
              }
            }
            w++;
        }
        if(str.back()=="")
          str.pop_back();                   //     for(int i=0; i<str.size(); i++) cout << str[i].c_str() << endl;
        //if(str.size() > 2)
        //{ //cout << str.size();
        //   cout << str[str.size() -1] << endl;
        //  printf("ERROR: string \" %s \": I can't recognize all the parameters\n", s.c_str());
        //  exit(1);
        //}
    
}



bool blank(string& s)
{
  size_t pos = s.length();
  size_t w = 0;
  bool blank = true;
  while(w!=pos) 
  {
      if(s[w]!= ' ')
      {
           blank = false;
           break;
      }
      w++;
  }
  return blank;
}

void End_Section(ifstream& fp_input, string str)
{
  string s;
  while((getline(fp_input,s)))
  {
   // cout << s << endl;

    if(0 == strcasecmp(s.c_str(), "")){}
    else
    {
          
          if(!blank(s))
          {
              size_t pos = s.find("}");
              if(pos == -1)
              { 
                  printf("Cannot find end of %s. Correct your input.\n", str.c_str());
                  exit(1);
              }
              else
              {
                 // printf("%s section read.\n", str.c_str());
                  return;
              } 
          }
    } 
  }
}




void read_laser_intensity(vector<string> & str, Laser& pulse){
  if(str.size()==2){ 
    double I = atof(str[1].c_str());
    pulse.set_IntensityWcm2(I); 
  }
  else if(0 == strcasecmp( str[2].c_str(), "wcm2") ){
    double I = atof(str[1].c_str());
    pulse.set_IntensityWcm2(I); 
  }
  else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(),  "atomicunit") ) {
    double I = atof(str[1].c_str())*intensity_au_Wcm2;
    pulse.set_IntensityWcm2(I); 
  }
  else {
    printf("Cannot understand %s as units.", str[2].c_str());
    exit(1);
  }
}







void read_laser_wavelenght(vector<string> & str, Laser& pulse){
  if(str.size()==2) {
    double lambda = atof(str[1].c_str());
    pulse.set_lambdanm(lambda); 
  }
  else if( 0 == strcasecmp( str[2].c_str(), "nm") ){
    double lambda = atof(str[1].c_str());
    pulse.set_lambdanm(lambda);
  }
  else if( 0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit") ) {      
    double lambda = atof(str[1].c_str())*10.*space_au_A;
    pulse.set_lambdanm(lambda);
  }
  else  {
    printf("Cannot understand %s as units.", str[2].c_str());
    exit(1);
  }
}






void set_laser_polarization(vector<string> & str, Laser& pulse){

  if(0 == strcasecmp( str[1].c_str(), "circular")){
    double value1 = atof(str[2].c_str()); //cout << value1 << endl;
    double value2 = atof(str[3].c_str()); //cout << value2 << endl;
    double value3 = atof(str[4].c_str()); //cout << value3 << endl;
    double value4 = atof(str[5].c_str()); //cout << value1 << endl;
    double value5 = atof(str[6].c_str()); //cout << value2 << endl;
    double value6 = atof(str[7].c_str()); //cout << value3 << endl;

    double value=sqrt(value1*value1 + value2*value2 + value3*value3); value1/=value; value2/=value; value3/=value;
    value=sqrt(value4*value4 + value5*value5 + value6*value6); value4/=value; value5/=value; value6/=value;

    pulse.set_pol(value1,value2,value3,value4,value5,value6);
  }
  else if(0 == strcasecmp( str[1].c_str(), "linear") || str.size() == 4 ) {
    double value1 = atof(str[1].c_str()); //cout << value1 << endl;
    double value2 = atof(str[2].c_str()); //cout << value2 << endl;
    double value3 = atof(str[3].c_str()); //cout << value3 << endl;

    double value=sqrt(value1*value1 + value2*value2 + value3*value3); value1/=value; value2/=value; value3/=value;
    double zero = 0.;
    pulse.set_pol(value1,value2,value3,zero,zero,zero);
  }
}





void set_laser_window(vector<string> & str, Laser& pulse, bool& pseudoPW){
  if(!pseudoPW) {
    printf("You defined a window for a pulse that is not pseudoPW!! Error.\n");
    exit(1);
  }

  double w1, w2;
  w1 = atof(str[1].c_str());
  w2 = atof(str[2].c_str());
  if(str.size() > 3 && 0 == strcasecmp( str[3].c_str(),  "fs")) {
    w1 *= time_fs_au;
    w2 *= time_fs_au;
  }

  pulse.set_window(w1,w2);

}







void set_laser_sigma(vector<string> & str, Laser& pulse, double & sigma, bool & sin2){
  if(sin2 == true){
    printf("You can't define sigma for sin2 wave. Change your input.\n");
    exit(1);
  }
  else {
    if(str.size()<3) {
      sigma = atof(str[1].c_str());
    }
    else if(0 == strcasecmp( str[2].c_str(), "fs")){//cout << "si " << endl;
      sigma = atof(str[1].c_str());
      sigma *= time_fs_au;
    }
    else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit")) {
      sigma = atof(str[1].c_str());
    }
    pulse.set_ncycle_sigma(sigma);
  }
}




void set_laser_ncycle(vector<string> & str, Laser& pulse, double & sigma, bool & gaussian){
  if(gaussian == true)  {
    printf("You can't define cycles for gaussian wave. Change your input.\n");
    exit(1);
  }

  sigma = atof(str[1].c_str());
  pulse.set_ncycle_sigma(sigma);
}


void read_laser_t0(vector<string> & str, double & t0){

    t0 = atof(str[1].c_str());
    if(0 == strcasecmp( str[2].c_str(), "fs"))  t0 *= time_fs_au;

}







void read_Laser_pump(vector<string> & str, ifstream& fp_input, Laser& pulse, 
  bool & gaussian, double & sigma, string & s, size_t & pos){

  bool pseudoPW = false;
  bool sin2 = false;
  double t0 = 0;

  gaussian = false;
  if(str.size() < 2){
    printf("in %s you have to specify laser type: gaussian or sin2. \n", s.c_str());
    exit(1);
  }

  if( 0 == strcasecmp(str[1].c_str(), "gaussian") )  gaussian = true;
  else if(str[1] == "sin2")  sin2 = true;
  else if(str[1] == "pw")  pseudoPW = true;
  else  {
    printf("Allowed types for laserpump: gaussian, sin2. Your type \" %s \" is not available.\n", str[1].c_str());
    exit(1);
  }

  bool envelope = false;
  pulse.set_boolean(envelope,gaussian,pseudoPW,sin2);

  if(!(getline(fp_input, s)))  {
    printf("Can't find the end of laserpump section.\n"); 
    exit(1);
  } 



  while(s==""){}

  while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), "")) {
    str.clear();
    pos = s.length();
    Separate_string(s, str, pos);
    if( 0 ==strcasecmp(str[0].c_str(), "frequency") ){//cout << str.size() << endl;
      if(str.size()==2){
        double wl = atof(str[1].c_str());
        pulse.set_freq(wl); //cout << s << endl;
      } 
      else {
        printf("Cannot understand %s as units.", str[2].c_str());
        exit(1);
      }
    } 
    else if ( 0 == strcasecmp( str[0].c_str(),  "phase") ){
      double phi=atof(str[1].c_str());
      pulse.set_phase(phi);
    }
    else if ( 0 == strcasecmp( str[0].c_str(),  "window") ){
      set_laser_window(str, pulse, pseudoPW);
    }
    else if ( 0 == strcasecmp( str[0].c_str(),  "wavelength") ){
      read_laser_wavelenght(str, pulse);
    }
    else if(0 == strcasecmp( str[0].c_str(),  "intensity") ) {
      read_laser_intensity(str, pulse);
    }
    else if(0 == strcasecmp( str[0].c_str(), "polarization") ) {  //cout << s <<" " << str.size() <<  endl;
      set_laser_polarization(str, pulse);
    }

    else if(0 == strcasecmp( str[0].c_str(), "sigma") ) {
      set_laser_sigma(str, pulse, sigma, sin2);
    }

    else if(0 == strcasecmp( str[0].c_str(), "cycles") ) {
      set_laser_ncycle(str, pulse,  sigma, gaussian);
    }
    else if(0 == strcasecmp( str[0].c_str(), "t0_pulse") ) {
      read_laser_t0(str, t0);
    }

    else if(!blank(str[0]) || (0 != strcasecmp(s.c_str(), ""))) {
      printf("cannot recognize key %s", s.c_str());
      exit(1);
    }

    if(!(getline(fp_input,s))) {
      printf("Can't find the end of laserpump section.\n"); 
      exit(1);
    }
  }//end of while

  if(gaussian) t0 += 4.*sigma;
  pulse.set_t0(t0);
  //    printf("Laser pump section read.\n");

}





void Read_Input
   (ifstream& fp_input, vec2d& a, vec2d& b, //unit cell
   vec1i& Nk, vec2d& kpt, bool& kptread,         //k space
    vec1i& Nb,   //bands
   //int& Nc, int& Nv, int& Nch,              //bands
   vector<Laser> & Laser_pumps, Laser& pulse2, double& DELAY,//Coord_B& u1, Coord_B& u2, bool& gaussian1, bool& gaussian2, double& sigma1, double& sigma2, double& DELAY,//lasers
   string& iMode, string& TBtype, double& dt, double& t_fin, //tdse options
   bool& iWFDs, int& wfd_resolution, Coulomb_parameters& Coulomb_set, bool& iCurrent, bool& iTAbs, bool& iTAbsK,    //observables
   double& T1, double& T2, double& Tch,
   int& nTAk, vec2d& TAkpt, vec1i& tagTAk,
    vec1d& hopping,
    double& FermiEnergy,
    methods_Diff_Eq &   Diff_Eq, int & it_resolution,   //decoherences
    double& change_gap_constant, bool & Vectorization)
{
  Coord_B u1; 
  Coord_B u2; 
  bool gaussian1; 
  bool gaussian2; 
  double sigma1; 
  double sigma2; 
  string s;
  double value;
  int i_pump = 0; 
  int num_pump_pulses = 1;
  Laser_pumps.resize(1);
  while(getline(fp_input,s))
  {
   // cout << s << endl;
   if(s!="")
   {
      vector<string> str;
      size_t pos = s.find("{");
    


      if (pos == -1)
      {
         printf("Can't find the character { in the string \" %s \". \nPlease, review your input.\n", s.c_str());
         exit(1);
      }


      else  Separate_string(s, str, pos);
      
      if ( 0 ==  strcasecmp(str[0].c_str(),"unitcell") )
      {
            if(!(fp_input >> a[0][0] >> a[0][1] >> a[0][2]
                     >> a[1][0] >> a[1][1] >> a[1][2]
                     >> a[2][0] >> a[2][1] >> a[2][2]))
              {printf("Cannot read unit cell because the file ended. Correct your input."); exit(1);}

            // cout << a[0][0] << a[0][1] << a[0][2]
            //         << a[1][0] << a[1][1] << a[1][2]
            //         << a[2][0] << a[2][1] << a[2][2];
            double value;
            VecProd(b[0],a[1],a[2]); VecProd(b[1],a[2],a[0]); VecProd(b[2],a[0],a[1]);
            value=DotProd(a[0],b[0])/2./pi; b[0][0]/=value; b[0][1]/=value; b[0][2]/=value;
            value=DotProd(a[1],b[1])/2./pi; b[1][0]/=value; b[1][1]/=value; b[1][2]/=value;
            value=DotProd(a[2],b[2])/2./pi; b[2][0]/=value; b[2][1]/=value; b[2][2]/=value;

            //cout << str.size() << endl;
            if (str.size() == 2) 
            {
              if( 0 == strcasecmp(str[1].c_str(), "angstrom"))
              {
                    for(int i=0; i<3; i++)
                    {
                      for(int j=0; j<3; j++)
                      {
                         a[i][j] *= space_A_au;
                         b[i][j] /= space_A_au;
                      } 
                    }
              }
              else if( 0==strcasecmp(str[1].c_str(), "au") || strcasecmp(str[1].c_str(), "atomicunit") ) {}
              else
              {
                 printf("Can't recognize the units %s in the line %s, please modify your input.\nAllowed units: angstrom, au / atomicunit", str[1].c_str(), s.c_str());
                 exit(1);
              }
            }
            Coord_B::set_crys_to_cart(b);

            End_Section(fp_input, str[0]);
 
      }




      else if(0 ==  strcasecmp( "nkpt", str[0].c_str() ) ) 
      {
          if(!(fp_input >> Nk[0] >> Nk[1] >> Nk[2]))                         {printf("Cannot read Nkpt because the file ended. Correct your input."); exit(1);}

          if (str.size() == 2)
          {
              if(0 ==  strcasecmp(str[1].c_str(), "read")) kptread = true;
              else printf("%s can't have units %s. Please modify your input.", str[0].c_str(), str[1].c_str());
          }
          if(kptread)
          {
              kpt.resize(Nk[0]*Nk[1]*Nk[2],3);
              for(int ik=0; ik<Nk[0]*Nk[1]*Nk[2];ik++)
              {
                  fp_input >> kpt[ik][0] >> kpt[ik][1] >> kpt[ik][2];
              }

          }
          
          End_Section(fp_input, str[0]);


      }



      else if(0 == strcasecmp(str[0].c_str(), "bands"))
      {
          s="";
          while(s==""){getline(fp_input, s);}
          
          while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
          {
              str.clear();
              pos = s.length();
              Separate_string(s, str, pos);
              if(strcasecmp(str[0].c_str(), "conduction")==0)
                Nb[2] = atoi(str[1].c_str()); //cout << "cond "<< str[1].c_str() << endl;}
              else if (strcasecmp(str[0].c_str(), "valence")==0)
                Nb[1] = atoi(str[1].c_str());   // cout << "valence "<< str[1].c_str() << endl;}
              else if (strcasecmp(str[0].c_str(), "core")==0)
                Nb[0] = atoi(str[1].c_str()); //cout << "core "<< str[1].c_str() << endl;}
              else if (strcasecmp(str[0].c_str(), "conduction_shift")==0)
                change_gap_constant = atof(str[1].c_str());
              else if(!blank(str[0]) || (0 != strcasecmp(str[0].c_str(), "")))
              {
                printf("Can't recognize band type %s. Possible band types: conduction, valence, core.\n", str[0].c_str());
                exit(1);
              }
              if(!(getline(fp_input,s))) 
              {
                  printf("Can't find the end of bands section.\n"); 
                  exit(1);
              }
          }
          
         /// printf("Bands section read.\n");
      }




      else if( 0 == strcasecmp(str[0].c_str(), "laserpump") and (i_pump < num_pump_pulses) ){
        Laser pulse_pump;
        read_Laser_pump(str, fp_input, pulse_pump, gaussian1, sigma1, s, pos);
        Laser_pumps[i_pump] = pulse_pump;// .push_back(pulse_pump);
        i_pump+= 1;
      }//end of laserpump

      else if( 0 == strcasecmp(str[0].c_str(), "laserprobe"))
      {
        gaussian2=false;
        bool sin2 = false;
        bool pseudoPW = false;
             if(str.size() < 2)
             {
                 printf("in %s you have to specify laser type: gaussian or sin2. \n", s.c_str());
                 exit(1);
             }
             bool envelope;
             if( 0 == strcasecmp(str[1].c_str(), "gaussian") )
             {
                gaussian2 = true;
              }
             else if(str[1] == "sin2")
             {
                  sin2 = true;
             }
             else if(str[1] == "pw")
             {
                pseudoPW = true;
             }
             else
             {
                printf("Allowed types for laserpump: gaussian, sin2. Your type \" %s \" is not available.\n", str[1].c_str());
                exit(1);
             } 
             
             if(!(getline(fp_input, s)))
             {
                  printf("Can't find the end of laserprobe section.\n"); 
                  exit(1);
             } 
             envelope=true;
             pulse2.set_boolean(envelope,gaussian2,pseudoPW,sin2);


             while(s==""){}
          
          while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
          {
              str.clear();
              pos = s.length();
              Separate_string(s, str, pos);
               if( 0 ==strcasecmp(str[0].c_str(), "frequency") )
               {//cout << str.size() << endl;
                 if(str.size()==2)
                 {
                   double wl = atof(str[1].c_str());
                   pulse2.set_freq(wl); //cout << s << endl;
                 }
                 else
                 {
                     printf("Cannot understand %s as units.", str[2].c_str());
                     exit(1);
                 }
               }              
               
               else if ( 0 == strcasecmp( str[0].c_str(),  "wavelength") )
               {
                    if(str.size()==2)
                    {
                        double lambda = atof(str[1].c_str());
                        pulse2.set_lambdanm(lambda); 
                    }
                    else if( 0 == strcasecmp( str[2].c_str(), "nm") )
                    {
                        double lambda = atof(str[1].c_str());
                        pulse2.set_lambdanm(lambda); 
                    }
                    else if( 0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit") )
                    {      
                        double lambda = atof(str[1].c_str())*10.*space_au_A;
                        pulse2.set_lambdanm(lambda); 
                    }
                    else
                    {
                          printf("Cannot understand %s as units.", str[2].c_str());
                          exit(1);
                    }
                }
               else if ( 0 == strcasecmp( str[0].c_str(),  "window") )
               {
                  if(!pseudoPW)
                  {
                    printf("You defined a window for a pulse that is not pseudoPW!! Error.\n");
                    exit(1);
                  }

                  double w1, w2;
                  w1 = atof(str[1].c_str());
                  w2 = atof(str[2].c_str());
                  if(str.size() > 3 && 0 == strcasecmp( str[3].c_str(),  "fs"))
                  {
                      w1 *= time_fs_au;
                      w2 *= time_fs_au;
                  }
                  pulse2.set_window(w1,w2);
               }
                else if(0 == strcasecmp( str[0].c_str(),  "intensity") )
                {
                    if(str.size()==2)
                    {
                        double I = atof(str[1].c_str());
                        pulse2.set_IntensityWcm2(I); 
                    }
                    else if(0 == strcasecmp( str[2].c_str(), "wcm2") )
                    {
                        double I = atof(str[1].c_str());
                        pulse2.set_IntensityWcm2(I); 
                    }
                    else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(),  "atomicunit") )
                    {
                        double I = atof(str[1].c_str())*intensity_au_Wcm2;
                        pulse2.set_IntensityWcm2(I); 
                    }
                    else
                    {
                          printf("Cannot understand %s as units.", str[2].c_str());
                          exit(1);
                    }
                }
                else if(0 == strcasecmp( str[0].c_str(), "polarization") )
                {  //cout << s <<" " << str.size() <<  endl;
                    if(0 == strcasecmp( str[1].c_str(), "circular"))
                    {
                        double value1 = atof(str[2].c_str()); //cout << value1 << endl;
                        double value2 = atof(str[3].c_str()); //cout << value2 << endl;
                        double value3 = atof(str[4].c_str()); //cout << value3 << endl;
                        double value4 = atof(str[5].c_str()); //cout << value1 << endl;
                        double value5 = atof(str[6].c_str()); //cout << value2 << endl;
                        double value6 = atof(str[7].c_str()); //cout << value3 << endl;

                        double value=sqrt(value1*value1 + value2*value2 + value3*value3); value1/=value; value2/=value; value3/=value;
                        value=sqrt(value4*value4 + value5*value5 + value6*value6); value4/=value; value5/=value; value6/=value;
                        
                        pulse2.set_pol(value1,value2,value3,value4,value5,value6);
                    }
                    else if(0 == strcasecmp( str[1].c_str(), "linear") || str.size() == 4 )
                    {
                        double value1 = atof(str[1].c_str()); //cout << value1 << endl;
                        double value2 = atof(str[2].c_str()); //cout << value2 << endl;
                        double value3 = atof(str[3].c_str()); //cout << value3 << endl;
  
                        double value=sqrt(value1*value1 + value2*value2 + value3*value3); value1/=value; value2/=value; value3/=value;
                        double zero = 0.;
                        pulse2.set_pol(value1,value2,value3,zero,zero,zero);
                    }
                }
                else if(0 == strcasecmp( str[0].c_str(), "sigma") )
                {
                    if(sin2 == true)
                    {
                        printf("You can't define sigma for sin2 wave. Change your input.\n");
                        exit(1);
                    }
                    else
                    {
                        if(str.size()<3)
                        {
                            sigma2 = atof(str[1].c_str());
                        }
                        else if(0 == strcasecmp( str[2].c_str(), "fs"))
                        {//cout << "si " << endl;
                            sigma2 = atof(str[1].c_str());
                            sigma2 *= time_fs_au;
                        }
                        else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit"))
                        {
                            sigma2 = atof(str[1].c_str());
                        }

                        pulse2.set_ncycle_sigma(sigma2);
                    }
                }
               



               else if( 0 ==strcasecmp(str[0].c_str(), "delay") )
               {//cout << str.size() << endl;
                 if(str.size()==2)
                 {
                   DELAY = atof(str[1].c_str())*time_fs_au; //cout << s << endl;
                 }
                 else if(0 ==strcasecmp(str[2].c_str(), "fs"))
                 {
                   DELAY = atof(str[1].c_str())*time_fs_au; //cout << s << endl;
                 }
                 else if(0 ==strcasecmp(str[2].c_str(), "au")||(0 ==strcasecmp(str[2].c_str(), "atomicunit")))
                 {
                   DELAY = atof(str[1].c_str()); //cout << s << endl;
                 }
                 else
                 {
                    printf("Error in reading time delay. Please review your input.\n");
                    exit(1);
                 }

               }              


                else if(0 == strcasecmp( str[0].c_str(), "cycles") )
                {
                    if(gaussian2 == true)
                    {
                        printf("You can't define cycles for gaussian wave. Change your input.\n");
                        exit(1);
                    }

                    sigma2 = atof(str[1].c_str());
                    pulse2.set_ncycle_sigma(sigma2);
                }
                else
                {
                    printf("cannot recognize key %s in laser probe.\n", s.c_str());
                }

              if(!(getline(fp_input,s))) 
              {
                  printf("Can't find the end of laserprobe section.\n"); 
                  exit(1);
              }

              double t0;


              }//end of while
              double t0;
              if(gaussian1) t0=4.*sigma1 + DELAY; 
              else t0=(Laser_pumps[0].tf - Laser_pumps[0].t0)/2. + DELAY - pulse2.ncycle/2.*pulse2.Period;
              pulse2.set_t0(t0);
        //      printf("Laser probe section read.\n");
          }//end of probe



          else if ( 0 == strcasecmp(str[0].c_str(), "TDSE"))
          {
             if(!(getline(fp_input, s)))
             {
                  printf("Can't find the end of tdse section.\n"); 
                  exit(1);
             } 


             while(s.find("}") == -1)
             {
                  str.clear();
                  pos =s.length();
                  Separate_string(s, str, pos);
                  if( 0 ==strcasecmp(str[0].c_str(), "Wannier") )
                  {
                      iMode = "W";
                  }//end of wannier
                  else if( 0 ==strcasecmp(str[0].c_str(), "Crystal") )
                  {
                      iMode = "CY";
                  }//end of crystal
                  else if (0 ==strcasecmp(str[0].c_str(), "tightbinding") )
                  {
                      if(str.size() < 2)
                      {
                          printf("Please specify the TightBinding type. Available types: Graphene, CoreGraphene, GrapheneZurron, BoronNitride, BN1b, DLGraphene, Vampa2015, Haldane_CoreBN_Bedge, gen2d_hexagonal, GeS, GeS_HSE06.\n");
                          exit(1);
                      }
                      else
                      {
                          TBtype = str[1].c_str();
                          if ( 0 == strcasecmp(TBtype.c_str(), "graphene") ||0 == strcasecmp(TBtype.c_str(), "GrapheneZurron") )
                          {
                            double LatticeConstant = 2.46*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(pi/6.);          a[1][2] = LatticeConstant*sin(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(pi/6.);          a[2][2] =-LatticeConstant*sin(pi/6.);

                              Nb[0] = 0;  Nb[1] = 1; Nb[2] = 1; 

                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "CoreGraphene") || 0 == strcasecmp(TBtype.c_str(), "Graphene_costD") || 0 == strcasecmp(TBtype.c_str(), "Graphene_symmetricD") )
                          {
                            double LatticeConstant = 2.46*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(pi/6.);          a[1][2] = LatticeConstant*sin(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(pi/6.);          a[2][2] =-LatticeConstant*sin(pi/6.);
                              
                               Nb[0] = 2;  Nb[1] = 1; Nb[2] = 1;
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "BoronNitride"))
                          {
                            double LatticeConstant = 2.5*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(pi/6.);          a[1][2] = LatticeConstant*sin(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(pi/6.);          a[2][2] =-LatticeConstant*sin(pi/6.);
                              
                              
                              Nb[0] = 0; Nb[1] = 1; Nb[2] = 1; 
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "CoreBN_Nedge") || 0 == strcasecmp(TBtype.c_str(), "CoreBN_Bedge") )
                          {
                            double LatticeConstant = 2.5*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(pi/6.);          a[1][2] = LatticeConstant*sin(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(pi/6.);          a[2][2] =-LatticeConstant*sin(pi/6.);
                              
                              
                               Nb[0] = 1;  Nb[1] = 1; Nb[2] = 1;
                          }

                          else if ( 0 == strcasecmp(TBtype.c_str(), "BN1b"))
                          {
                            double LatticeConstant = 2.5*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(pi/6.);          a[1][2] = LatticeConstant*cos(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(pi/6.);          a[2][2] =-LatticeConstant*cos(pi/6.);
                              
                              Nb[0] = 0; Nb[1] = 1; Nb[2] = 0;
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "DLGraphene"))
                          {
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;          a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = 4.0257407400;          a[1][2] =-2.3242630000;
                              a[2][0] = 0.0000000000;       a[2][1] = 4.0257407400;          a[2][2] = 2.3242630000;
                              
                              Nb[0] = 0; Nb[1] = 2; Nb[2] = 2;
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "Vampa2015"))
                          {
                              a[0][0] = 5.3200000000;       a[0][1] = 0.0000000000;          a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = 6.1400000000;          a[1][2] = 0.0000000000;
                              a[2][0] = 0.0000000000;       a[2][1] = 0.0000000000;          a[2][2] = 9.8300000000;
                              
                              Nb[0] = 0; Nb[1] = 1; Nb[2] = 1;
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "GeS") || 0 == strcasecmp(TBtype.c_str(), "GeS_HSE06") )
                          {
                              double ax = (4.53/2.)*space_A_au;
                              double ay = (3.63/2.)*space_A_au;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = ax;          a[1][2] = ay;
                              a[2][0] = 0.0000000000;       a[2][1] = ax;          a[2][2] =-ay;
                              
                              Nb[0] = 0; Nb[1] = 1; Nb[2] = 1;
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "Haldane_CoreBN_Bedge") || 0 == strcasecmp(TBtype.c_str(), "Haldane_CoreBN_Nedge") )
                          {
                            double LatticeConstant = 2.5*space_A_au;
                            double Detune_angle = 1.0;
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                        a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = LatticeConstant*cos(Detune_angle * pi/6.);          a[1][2] = LatticeConstant*sin(Detune_angle * pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = LatticeConstant*cos(Detune_angle * pi/6.);          a[2][2] =-LatticeConstant*sin(Detune_angle * pi/6.);
                              
                              
                              Nb[0] = 1; Nb[1] = 1; Nb[2] = 1;
                              
                              hopping.resize(4);
                              
                              for (int i=2; i<str.size(); i++) {
                                  if(0 ==strcasecmp(str[i].c_str(), "gap="))
                                      hopping[0] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "t1="))
                                      hopping[1] = atof(str[i+1].c_str())*energy_eV_au;
                                  if(0 ==strcasecmp(str[i].c_str(), "t2="))
                                      hopping[2] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "phi0="))
                                      hopping[3] = atof(str[i+1].c_str());
                              }
                          }
                          else if ( 0 == strcasecmp(TBtype.c_str(), "gen2d_hexagonal") )
                          {
                              hopping.resize(8);
                              for (int i=2; i<str.size(); i++) {
                                  if(0 ==strcasecmp(str[i].c_str(), "a="))
                                      hopping[0] = atof(str[i+1].c_str())*space_A_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "gap="))
                                      hopping[1] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "t1="))
                                      hopping[2] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "t2="))
                                      hopping[3] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "phi0="))
                                      hopping[4] = atof(str[i+1].c_str());
                                  else if(0 ==strcasecmp(str[i].c_str(), "core="))
                                      hopping[5] = atof(str[i+1].c_str())*energy_eV_au;
                                  else if(0 ==strcasecmp(str[i].c_str(), "chv_dipole="))
                                      hopping[6] = atof(str[i+1].c_str());
                                  else if(0 ==strcasecmp(str[i].c_str(), "chc_dipole="))
                                      hopping[7] = atof(str[i+1].c_str());
                              }
                              a[0][0] = 6.3303100000;       a[0][1] = 0.0000000000;                   a[0][2] = 0.0000000000;
                              a[1][0] = 0.0000000000;       a[1][1] = hopping[0]*cos(pi/6.);          a[1][2] = hopping[0]*sin(pi/6.);
                              a[2][0] = 0.0000000000;       a[2][1] = hopping[0]*cos(pi/6.);          a[2][2] = -hopping[0]*sin(pi/6.);
                              
                              
                              Nb[0] = 1; Nb[1] = 1; Nb[2] = 1;
                          }


                          
                      }
                      double value;
                      VecProd(b[0],a[1],a[2]); VecProd(b[1],a[2],a[0]); VecProd(b[2],a[0],a[1]);
                      value=DotProd(a[0],b[0])/2./pi; b[0][0]/=value; b[0][1]/=value; b[0][2]/=value;
                      value=DotProd(a[1],b[1])/2./pi; b[1][0]/=value; b[1][1]/=value; b[1][2]/=value;
                      value=DotProd(a[2],b[2])/2./pi; b[2][0]/=value; b[2][1]/=value; b[2][2]/=value;
                      Coord_B::set_crys_to_cart(b);

                      iMode = "TB";
                  }//end of tightbinding
                  else if (0 ==strcasecmp(str[0].c_str(), "FermiEnergy") )
                  {
                      if(str.size() < 3)
                          FermiEnergy = atof(str[1].c_str());
                      else
                          exit(1);

                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "dt") )
                  {
                      if(str.size() < 2)
                      {
                          dt = atof(str[1].c_str());
                      }
                      else if(0 ==strcasecmp(str[2].c_str(), "fs"))
                      {
                          dt = atof(str[1].c_str())*time_fs_au;
                      }
                      else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[0].c_str(), "atomicunit")))
                      {
                         dt = atof(str[1].c_str());
                      }                         
                      else
                      {
                          printf("Units %s aren't allowed for dt.\n", str[2].c_str());
                          exit(1);
                      }
                   }
                  else if(0 ==strcasecmp(str[0].c_str(), "t_fin") )
                  {
                      if(str.size() < 2)
                      {
                          t_fin = atof(str[1].c_str());
                      }
                      else if(0 ==strcasecmp(str[2].c_str(), "fs"))
                      {
                          t_fin = atof(str[1].c_str())*time_fs_au;
                      }
                      else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[0].c_str(), "atomicunit")))
                      {
                         t_fin = atof(str[1].c_str());
                      }                         
                      else
                      {
                          printf("Units %s aren't allowed for dt.\n", str[2].c_str());
                          exit(1);
                      }
                   }
                  else if(0 ==strcasecmp(str[0].c_str(), "epsStepAbs") )
                  {
                    Diff_Eq.epsStepAbs = atof(str[1].c_str());
         
                   }
                  else if(0 ==strcasecmp(str[0].c_str(), "dynamical_dt_evolution") ) {
                    Diff_Eq.const_dt_evolution = false;
                    Diff_Eq.dynamical_dt_evolution = true; 
                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "No_Vectorization") ) {
                    Vectorization = false;
                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "Taylor") ) {
                    Diff_Eq.Taylor= true; 
                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "num_pump_pulses") ) {
                    num_pump_pulses = atof(str[1].c_str());
                    Laser_pumps.resize(num_pump_pulses);
                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "TaylorOrder") ) {
                    Diff_Eq.TaylorOrder= atof(str[1].c_str()); 
                  }
                  else if(0 ==strcasecmp(str[0].c_str(), "Gap_correction") ) {
                    Diff_Eq.Gap_correction= atof(str[1].c_str()); 
                  }


                  if(!(getline(fp_input, s)))
                  {
                       printf("Can't find the end of tdse section.\n"); 
                       exit(1);
                  } 


             }//end of reading tdse
          //   printf("TDSE section read.\n");

          }//end of tdse


          else if ( 0 == strcasecmp(str[0].c_str(), "observables"))
          {
             if(!(getline(fp_input, s)))
             {
                  printf("Can't find the end of observables section.\n"); 
                  exit(1);
             } 


             while(s.find("}") == -1)
             {
                  pos =s.length();
                  str.clear();
                  Separate_string(s, str, pos);
                  if(!blank(str[0]) && (0 != strcasecmp(str[0].c_str(), "")))
                  {
                    if      ( 0 ==strcasecmp(str[0].c_str(), "wfd") ) 
                    {
                        iWFDs = true;
                        if(str.size()==1 || ( 0 ==strcasecmp(str[0].c_str(), "1") )) wfd_resolution = 1;
                				else if(str.size() == 2) wfd_resolution = atof(str[1].c_str());
                    }
                    else if ( 0 ==strcasecmp(str[0].c_str(), "Current") ) iCurrent = true;
                    else if ( 0 ==strcasecmp(str[0].c_str(), "it_resolution") ) it_resolution = atof(str[1].c_str()); 
                    else if ( 0 ==strcasecmp(str[0].c_str(), "PrintPopulation") ){
                      Diff_Eq.PrintPopulation = true;
                      if (str.size() > 1){
                        Diff_Eq.start_print_time = atof(str[1].c_str()); // In a.u.!
                        Diff_Eq.end_print_time = atof(str[2].c_str());
                        Diff_Eq.step_print = atof(str[3].c_str());
                      }
                    } 
                    else if ( 0 ==strcasecmp(str[0].c_str(), "Coulomb_band_reconstruction") )   Coulomb_set.Coulomb_calc = true;
                    else if ( 0 ==strcasecmp(str[0].c_str(), "TransientAbsorption") || 0 ==strcasecmp(str[0].c_str(), "TAbs") ) iTAbs = true;
                    else if ( 0 ==strcasecmp(str[0].c_str(), "kresolvedTAbs") || 0 ==strcasecmp(str[0].c_str(), "kTAbs") || 0 ==strcasecmp(str[0].c_str(), "kTransientAbsorption") ) iTAbsK = true;
                    else
                    {
                        printf("option %s not available in the observables.\n", str[0].c_str());
                        exit(1);
                    }
                  } 
                  if(!(getline(fp_input, s)))
                  { 
                      printf("Can't find the end of observables section.\n"); 
                      exit(1);
                  }

             }
           }//end of observables
           else if ( 0 == strcasecmp(str[0].c_str(), "TAbsK"))
           {
               if(!iTAbsK)
               {
                  printf("Error in the input. You defined the section TAbsK without putting it in the observables.\n");
                  printf("Maybe you printed this section before the observables?\n");
                  exit(1);
               }
               s="";
               while(s==""){getline(fp_input, s);}
               str.clear();   //cout << s << endl;
               pos = s.length();
               Separate_string(s, str, pos);
               if(str.size() != 1)
               {
                  printf("Error in TAbsK. you didn't define in the first line the number of k points.");
                  exit(1);
               }
               else nTAk = atoi(str[0].c_str());
               TAkpt.resize(nTAk,3);
               tagTAk.resize(nTAk);

               for(int i=0; i<nTAk; i++)    
               {  
                  getline(fp_input, s); 
                  str.clear();
                  pos = s.length();
                  Separate_string(s, str, pos);
                  if(str.size() != 3)
                  {
                    printf("Error in defining k points in TAbsK. Please review your input.\n");
                    exit(1);
                  } 
                  for(int j=0; j<3; j++) TAkpt[i][j] = atof(str[j].c_str());
               }
               End_Section(fp_input, "TAbsK");
           }// end TAbsK
                 ////////////////////////////////////////////////////////////
          else if ( 0 == strcasecmp(str[0].c_str(), "Coulomb"))
          {
            Coulomb_set.Coulomb_calc = true;
            if(!(getline(fp_input, s))){
                 printf("Can't find the end of Coulomb section.\n");
                 exit(1);
            }
            while(s.find("}") == -1){
              str.clear();
              pos =s.length();
              Separate_string(s, str, pos);

              if ( 0 ==strcasecmp(str[0].c_str(), "qTF") ){
                Coulomb_set.qTF = atof(str[1].c_str());
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "r0") ){ // ONLY a.u.!
                Coulomb_set.r0 = atof(str[1].c_str());
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "G_distance") ){ // is size of BZ!
                Coulomb_set.G_distance = atof(str[1].c_str());
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "Rytova_Keldysh") ){
                Coulomb_set.Rytova_Keldysh = true;
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "epsilon_static") ){
                Coulomb_set.epsilon_static = atof(str[1].c_str());
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "Coulomb_diag_basis") ){
                Coulomb_set.Wannie_basis = false;
                Coulomb_set.Diagonal_basis = true;
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "Read_Coulomb_from_files") ){
                Coulomb_set.Calculate = false; // By default we calculate caoefficients, but we can also read from file
                Coulomb_set.Read_from_files = true;
              }

              if ( 0 ==strcasecmp(str[0].c_str(), "labelInput") ){
               Coulomb_set.labelInput = str[1].c_str();
              }
              
              if ( 0 ==strcasecmp(str[0].c_str(), "Print_new_band_dispersion") ){
               Coulomb_set.Print_new_band_dispersion = true;
              }
              if ( 0 ==strcasecmp(str[0].c_str(), "Ncut") ){
                Coulomb_set.Ncut = atof(str[1].c_str());
              }

              if(!(getline(fp_input, s)))                 {
                   printf("Can't find the end of Coulomb section.\n");
                   exit(1);
              }
            }

          } // end of Coulomb reading


          else if ( 0 == strcasecmp(str[0].c_str(), "decoherence"))
          {
             if(!(getline(fp_input, s)))
             {
                  printf("Can't find the end of decoherence section.\n"); 
                  exit(1);
             } 
             

             while(s.find("}") == -1)
             {
                      str.clear();
                      pos =s.length();
                      
                      Separate_string(s, str, pos);
                      if(str.size()<2)
                      {
                          printf("You have to specify a number for the decoherence term.");
                          exit(1);
                      }
                      if      ( 0 ==strcasecmp(str[0].c_str(), "diagonal") )
                      {
                         T1 = atof(str[1].c_str());
                         if(str.size() == 3)
                         {
                            if( 0 == strcasecmp(str[2].c_str(), "meV"))
                            {
                                T1*= energy_eV_au/1000.;
                            }
                            else if(0 == strcasecmp(str[2].c_str(), "au") || 0 == strcasecmp(str[2].c_str(), "atomicunit"))
                            {}
                            else
                            {
                                printf("Can't understand units %s in decoherence.\n", str[2].c_str());
                                exit(1);
                            } 
                          }
                        }
                        else if ( 0 ==strcasecmp(str[0].c_str(), "offdiagonal") )
                        {
                            T2 = atof(str[1].c_str());
                            if(str.size() == 3)
                            {
                                  if( 0 == strcasecmp(str[2].c_str(), "meV"))
                                  {
                                      T2*= energy_eV_au/1000.;
                                  }
                                  else if(0 == strcasecmp(str[2].c_str(), "au") || 0 == strcasecmp(str[2].c_str(), "atomicunit"))
                                  {}
                                  else
                                  {
                                    printf("Can't understand units %s in decoherence.\n", str[2].c_str());
                                    exit(1);
                                  } 
                            }
                        }
                        else if ( 0 ==strcasecmp(str[0].c_str(), "corehole") )    
                        {
                            Tch= atof(str[1].c_str());

                            if(str.size() == 3)
                            {
                              if( 0 == strcasecmp(str[2].c_str(), "meV"))
                              {
                                  Tch*= energy_eV_au/1000.;
                              }
                              else if(0 == strcasecmp(str[2].c_str(), "au") || 0 == strcasecmp(str[2].c_str(), "atomicunit"))
                              {}
                              else
                              {
                                  printf("Can't understand units %s in decoherence.\n", str[2].c_str());
                                  exit(1);
                              }
                            } 
                        }
                        else
                        {
                          printf("Unknow key %s in decoherence.", str[0].c_str()); 
                          exit(1);
                        }
                      
                      if(!(getline(fp_input, s)))
                      {
                           printf("Can't find the end of tdse section.\n"); 
                           exit(1);
                      } 
                    } 
          }//end of decoherence


   }//end of if       
 }//end of while

}



