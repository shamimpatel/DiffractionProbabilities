#include <iostream>
#include "Vector.h"
#include "math.h"
#include "FormFactorData.h"
#include "LatticePlane.h"
#include "AbsorbCoeffData.h"
#include <map>
#include "FileReading.h"
#include "ccd.h"
using namespace std;


float UpperLimit( float Value, float Limit)
{
    if( Value > Limit )
    {
        return Limit;
    }
    else
    {
        return Value;
    }
}

bool ApproxEqual( float a, float b, float tolerance)
{
    return fabs( a - b ) <= tolerance;
}

struct hklSet
{
    int h,k,l;
    int M; //Multiplicity (convenient to just have it here)
    hklSet( int h, int k, int l, int M)
    {
        this->h = h;
        this->k = k;
        this->l = l;
        this->M = M;
    }
};


//this is terrible programming practice....
const int _G_ISBCC = 0; //set to zero for BCC, 1 for FCC

bool isSetInList( hklSet Set, std::vector<hklSet> *List)
{
    for( int i=0; i < int(List->size()); i++ )
    {
        if( Set.h == ((*List)[i]).h && Set.k == ((*List)[i]).k && Set.l == ((*List)[i]).l)
        {
            return true;
        }
    }

    if( (Set.h + Set.k + Set.l)%2 != _G_ISBCC) //if not BCC/FCC, say it's in the list (though it isn't)
    {
        return true;
    }

    return false;
}

int main ()
{
    
    ifstream datafile("InputScript.txt");
    if(datafile.is_open() == false)
    {
        cout << "Error: Failed to open InputScript.txt" << endl;
        exit(1);
    }
    std::map<std::string,std::string> InputData;
    AddToMapFromFile(datafile, InputData);
    datafile.close();
////////////////////////////////////////////////////////////////////////////////
//          Generate Diffraction Probabilities
////////////////////////////////////////////////////////////////////////////////
    
    FormFactorData FormFactor( 0.0f, 2.0f, 10000);
    std::string FormFacFilename;
    StringFromMap("FormFactorData", InputData, FormFacFilename);
    //FormFactor.LoadData("TaFormFactor.txt");
    FormFactor.LoadData(FormFacFilename.c_str());

    AbsorbCoeffData MuData( 1.0f, 15.0f, 2000);
    std::string MuDataFilename;
    StringFromMap("AbsorptionData", InputData, MuDataFilename);
    //MuData.LoadData("TaAbsorbCoeff.txt");
    MuData.LoadData(MuDataFilename.c_str());

    float a0 = 3.31f; //lattice constant 2.87 angstroms (iron) 3.31 (Ta)

    Vector a1(1,0,0);
    Vector a2(0,1,0);
    Vector a3(0,0,1);

    a1 = a1*a0;
    a2 = a2*a0;
    a3 = a3*a0;

    float UnitCellVol = a2.Cross(a3).Dot(a1);

    Vector b1 = (a2.Cross(a3))/UnitCellVol;
    Vector b2 = (a3.Cross(a1))/UnitCellVol;
    Vector b3 = (a1.Cross(a2))/UnitCellVol;

    cout << "UnitCellVol:\t" << UnitCellVol << endl;

    

    ofstream BraggAngleFile;
    BraggAngleFile.open( "BraggAngle.txt" );
    if(BraggAngleFile.is_open() == false)
    {
        cout << "Failed to open BraggAngle.txt" << endl;
        exit(1);
    }

    ofstream ScatterProbFile;
    ScatterProbFile.open( "ScatterProb.txt" );
    if(ScatterProbFile.is_open() == false)
    {
        cout << "Failed to open ScatterProb.txt" << endl;
        exit(1);
    }

    ofstream RockingCurveFile;
    RockingCurveFile.open( "RockingCurve.txt" );
    if(RockingCurveFile.is_open() == false)
    {
        cout << "Failed to open RockingCurve.txt" << endl;
        exit(1);
    }

    ScatterProbFile << "Energy\t";
    BraggAngleFile  << "Energy\t";
    RockingCurveFile << "Energy\t";

    int indexMax = 5;

    std::vector< hklSet > hklPlanes;

    hklSet S(0,0,0,0);

    for( int n1 = 0; n1 <= indexMax; n1++)
    {
        if(n1 != 0)
        {
            S = hklSet(0 ,0 ,n1, 6);
            if( !isSetInList(S, &hklPlanes) )
            {
                hklPlanes.push_back(S);
                ScatterProbFile << 0  << 0  << n1 << "\t";
                BraggAngleFile << 0  << 0  << n1 << "\t";
                RockingCurveFile << 0  << 0  << n1 << "\t";
            }
            S = hklSet(0 ,n1,n1,12);
            if( !isSetInList(S, &hklPlanes) )
            {
                hklPlanes.push_back(S);
                ScatterProbFile << 0  << n1 << n1 << "\t";
                BraggAngleFile << 0  << n1 << n1 << "\t";
                RockingCurveFile << 0  << n1 << n1 << "\t";
            }
            S = hklSet(n1,n1,n1, 8);
            if( !isSetInList(S, &hklPlanes) )
            {
                hklPlanes.push_back(S);
                ScatterProbFile << n1 << n1 << n1 << "\t";
                BraggAngleFile << n1 << n1 << n1 << "\t";
                RockingCurveFile << n1 << n1 << n1 << "\t";
            }
        }

        for( int n2 = n1; n2 <= indexMax; n2++)
        {
            if( n1 != 0 && n2 !=0)
            {
                S = hklSet(0 ,n1,n2,24);
                if( !isSetInList(S, &hklPlanes) )
                {
                    hklPlanes.push_back(S);
                    ScatterProbFile << 0  << n1 << n2 << "\t";
                    BraggAngleFile << 0  << n1 << n2 << "\t";
                    RockingCurveFile << 0  << n1 << n2 << "\t";
                }
                S = hklSet(n1,n1,n2,24);
                if( !isSetInList(S, &hklPlanes) )
                {
                    hklPlanes.push_back(S);
                    ScatterProbFile << n1 << n1 << n2 << "\t";
                    BraggAngleFile << n1 << n1 << n2 << "\t";
                    RockingCurveFile << n1 << n1 << n2 << "\t";
                }
            }

            for( int n3 = n2; n3 <= indexMax; n3++)
            {
                if( n1 != 0 && n2 !=0 && n3 != 0)
                {
                    S = hklSet(n1,n2,n3,48);
                    if( !isSetInList(S, &hklPlanes) )
                    {
                        hklPlanes.push_back(S);
                        ScatterProbFile << n1 << n2 << n3 << "\t";
                        BraggAngleFile << n1 << n2 << n3 << "\t";
                        RockingCurveFile << n1 << n2 << n3 << "\t";
                    }
                }
            }
        }
    }

    //ScatterProbFile << "Sum_" << indexMax << endl;
    ScatterProbFile << endl;
    BraggAngleFile << endl;
    RockingCurveFile << endl;

    
    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
    double DeltaE = 0.0005;//0.0005f;
    
    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    DoubleFromMap("DeltaE", InputData, DeltaE);
    
    int nEPoints = (MaxE - MinE)/DeltaE;

    //for( float Energy = 3.0f; Energy <= 10.0f; Energy += 0.001f)
    for( int EnergyTick = 0; EnergyTick < nEPoints; EnergyTick++)
    {
        double Energy = MinE + DeltaE*EnergyTick;

        ScatterProbFile << Energy << "\t";
        BraggAngleFile << Energy << "\t";
        RockingCurveFile << Energy << "\t";

        float Sum = 0.0f;

        for( int i=0; i<int(hklPlanes.size()); i++)
        {
            LatticePlane Plane( b1, b2, b3, hklPlanes[i].h,  hklPlanes[i].k,  hklPlanes[i].l, &FormFactor, UnitCellVol,
                               &MuData, hklPlanes[i].M);
            BraggAngleFile << Plane.FindBraggReflectionAngle( EnergyToWavelength(Energy) ) << "\t";
            float I = Plane.CalculatePowderScatter( EnergyToWavelength(Energy) );
            Sum += I;
            ScatterProbFile << I << "\t";
            RockingCurveFile << Plane.GaussianScherrerWidth(EnergyToWavelength(Energy), 1000) << "\t"; //1000A for grain size (physical?)
        }

        //ScatterProbFile << Sum;
        ScatterProbFile << endl;
        BraggAngleFile << endl;
        RockingCurveFile << endl;
    }

    BraggAngleFile.close();
    ScatterProbFile.close();
    RockingCurveFile.close();
    
////////////////////////////////////////////////////////////////////////////////
//          Comput Min/Max Bragg Angles
////////////////////////////////////////////////////////////////////////////////
    
  
    CCD CCDCamera = GenerateCCDFromInputScript("InputScript.txt");
    
    Vector CCDCorners[4];
    
    CCDCamera.GenerateCCDCorners( CCDCorners[0], CCDCorners[1], CCDCorners[2], CCDCorners[3]);
    
    double CCDXMin = CCDCorners[0].x, CCDXMax = CCDCorners[0].x, CCDYMin = CCDCorners[0].y, CCDYMax = CCDCorners[0].y;
    
    cout << "CCDBounds:" << endl;    
    for(int i = 0; i<4; i++)
    {
        cout << i << ":\t"; CCDCorners[i].Print();
        if(CCDCorners[i].x < CCDXMin)
        {
            CCDXMin = CCDCorners[i].x;
        }
        if(CCDCorners[i].y < CCDYMin)
        {
            CCDYMin = CCDCorners[i].y;
        }
        if(CCDCorners[i].x > CCDXMax)
        {
            CCDXMax = CCDCorners[i].x;
        }
        if(CCDCorners[i].y > CCDYMax)
        {
            CCDYMax = CCDCorners[i].y;
        }
        
    }
    
    Vector CrystalOrigin(0,-0.1,0);
    double CrystalXLength = 0.1;
    double CrystalYLength = 0.2;
    
    VectorFromMap("CrystalOrigin",InputData, CrystalOrigin);
    DoubleFromMap("CrystalXLength", InputData, CrystalXLength);
    DoubleFromMap("CrystalYLength", InputData, CrystalYLength);
    
    cout << "CrystalOrigin:\t"; CrystalOrigin.Print();
    cout << "CrystalDimensions:\tX:\t" << CrystalXLength << "\tY:\t" << CrystalYLength << endl;
    
    Vector CrystalCorners[4];
    
    CrystalCorners[0] = Vector(CrystalOrigin);
    CrystalCorners[1] = Vector(0,CrystalOrigin.y+CrystalYLength,0);
    CrystalCorners[2] = Vector(CrystalOrigin.x+CrystalXLength,CrystalOrigin.y,0);
    CrystalCorners[3] = Vector(CrystalOrigin.x+CrystalXLength,CrystalOrigin.y+CrystalYLength,0);
    
    cout << "Crystal Corners:" << endl;
    for(int i = 0; i<4; i++)
    {
        CrystalCorners[i].Print();
    }
    
    Vector Source( -3.0f, 0.0f, 5.0);    
    VectorFromMap("Source", InputData, Source);
    
    double MinBraggAngle = 10;
    double MaxBraggAngle = -1;
    
    //TODO: loop across crystal edges AND CCD edges??
    //adding in a 2% error should handle that? (tolerance now taken via input script)
    
    for(int iCrystalCorner = 0; iCrystalCorner < 4; iCrystalCorner++)
    {
        Vector SourceToCrystal = CrystalCorners[iCrystalCorner] - Source;
        SourceToCrystal = SourceToCrystal.Normalized();
        for(int iCCDCorner = 0; iCCDCorner < 4; iCCDCorner ++)
        {
            Vector CrystalToCCD = CCDCorners[iCCDCorner] - CrystalCorners[iCrystalCorner];
            CrystalToCCD = CrystalToCCD.Normalized();
            double BraggAngle = 0.5*acos( CrystalToCCD.Dot(SourceToCrystal));//0.5 as the calculation gives scattering angle;
            if( BraggAngle < MinBraggAngle )
            {
                MinBraggAngle = BraggAngle;
            }
            if( BraggAngle > MaxBraggAngle )
            {
                MaxBraggAngle = BraggAngle;
            }
        }
    }
    
    if(MinBraggAngle > 0.5*PI || MinBraggAngle < 0.0||
       MaxBraggAngle > 0.5*PI || MaxBraggAngle < 0.0)
    {
        cout << "Error: Calculated Bragg angles out of bounds!" << endl;
        exit(1);
    }
    
    double BraggAngleTolerance = 0.01;    
    DoubleFromMap("BraggAngleTolerance", InputData, BraggAngleTolerance);
    
    //This may push the angle limits out of [0,pi/2] but we don't really care if this happens since they're used to
    //optimise the scattering probability data. It is checked above since there must be something wrong for the
    //calculation to produce something out of the [0,pi/2] interval.
    MinBraggAngle -= MinBraggAngle*BraggAngleTolerance;
    MaxBraggAngle += MaxBraggAngle*BraggAngleTolerance;
    
    cout << "MinAngle:\t" << MinBraggAngle << "\tMaxAngle:\t" << MaxBraggAngle << endl;
    
    ofstream BraggAngleLimitsFile;
    BraggAngleLimitsFile.open( "BraggAngleLimits.txt" );
    
    //keep the first part as one word!
    BraggAngleLimitsFile << "MinimumBraggAngle:\t" << MinBraggAngle << endl;
    BraggAngleLimitsFile << "MaximumBraggAngle:\t" << MaxBraggAngle << endl;
    
    BraggAngleLimitsFile.close();
    
    return 0;

}
