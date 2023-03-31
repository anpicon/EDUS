//loop for printing matrices
ofstream pU, pH, pE, pDx, pDy, pDz;
pU.open("Output/U"  +srank.str()+".txt");
pH.open("Output/H"  +srank.str()+".txt");
pE.open("Output/E"  +srank.str()+".txt");
pDx.open("Output/Dx"+srank.str()+".txt");
pDy.open("Output/Dy"+srank.str()+".txt");
pDz.open("Output/Dz"+srank.str()+".txt");

for(int ikx=-100; ikx<100; ikx++)
{
    double kx = 4.*double(ikx)/100.;
    for(int iky=-100; iky<100; iky++)
    {
        double ky=4.*double(iky)/100.;
        kprint.setcart(kx,ky,0.);
        WModel.energy_U(Hk,Uk,kprint);
        Ek.fill(0.);
        pU << kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  "; 
        pH << kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  ";
        pDx<< kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  ";
        pDy<< kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  ";
        pDz<< kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  ";
        pE<< kprint.cart[0] << " " << kprint.cart[1] << " " << kprint.cart[2] << "  ";
        for(int ic=0; ic<Ncv; ic++)
        {
            for(int jc=0; jc<Ncv; jc++)
            {   
                pU  << real(Uk[ic][jc]) << " " << imag(Uk[ic][jc]) << " ";
                pH  << real(Hk[ic][jc]) << " " << imag(Hk[ic][jc]) << " ";
                pDx << real(WModel.dipolex(ic,jc,kprint)) << " " << imag(WModel.dipolex(ic,jc,kprint)) << " ";
                pDy << real(WModel.dipoley(ic,jc,kprint)) << " " << imag(WModel.dipoley(ic,jc,kprint)) << " ";
                pDz << real(WModel.dipolez(ic,jc,kprint)) << " " << imag(WModel.dipolez(ic,jc,kprint)) << " ";
                for(int kc=0; kc<Ncv; kc++)
                    for(int lc=0; lc<Ncv; lc++)
                        Ek[ic][jc] += Uk[ic][kc]*Hk[kc][lc]*conj(Uk[jc][lc]);
                pE << real(Ek[ic][jc]) << " " << imag(Ek[ic][jc]) << " ";
            }
        }
        pU  << endl;
        pH  << endl;
        pDx << endl;
        pDy << endl;
        pDz << endl;
        pE << endl;
    }
}

pU.close();
pH.close();
pDx.close();
pDy.close();
pDz.close();
pE.close();
