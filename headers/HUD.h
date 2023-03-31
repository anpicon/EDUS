TightBinding TB_Model;
Wannier WModel;
Crystal CYModel;

void gHUD(vec2x& H,vec2x& U,vec3x& D, Coord_B& k,string& iMode)
{
    if(iMode=="CY"){
        CYModel.energy_U(H,U,k);
        CYModel.dipole(D,k);
    }
    else if(iMode=="W"){
        WModel.energy_U(H,U,k);
        WModel.dipole(D,k);
    }
    else {
        TB_Model.energy_U(H,U,k);
        TB_Model.dipole(D,k);
    }
}

void gHU(vec2x& H,vec2x& U, Coord_B& k,string& iMode)
{
    if(iMode=="CY"){
        CYModel.energy_U(H,U,k);
    }
    else if(iMode=="W"){
        WModel.energy_U(H,U,k);
    }
    else {
        TB_Model.energy_U(H,U,k);
    }
}


