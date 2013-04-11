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

//this is terrible programming practice....
const int _G_ISBCC = 0; //set to zero for BCC, 1 for FCC

bool isPlaneInList( LatticePlane Plane, std::vector<LatticePlane> *List)
{
    for( int i=0; i < int(List->size()); i++ )
    {
        if( Plane.h == ((*List)[i]).h && Plane.k == ((*List)[i]).k && Plane.l == ((*List)[i]).l)
        {
            return true;
        }
    }

    if( (Plane.h + Plane.k + Plane.l)%2 != _G_ISBCC) //if not BCC/FCC, say it's in the list (though it isn't)
    {
        return true;
    }

    return false;
}

int main ()
{
    cout << endl;
    cout << "//////////////////////////////////////////////////////" << endl;
    cout << "//Now running Diffraction Probabiities Pre-Processor//" << endl;
    cout << "//////////////////////////////////////////////////////" << endl;
    cout << endl;
        
    
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
    
    
    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
    double DeltaE = 0.0005;//0.0005f;
    
    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    DoubleFromMap("DeltaE", InputData, DeltaE);
    
    int nEPoints = (MaxE - MinE)/DeltaE;    
    
    
    /*std::string FormFacFilename;
    StringFromMap("FormFactorData", InputData, FormFacFilename);
    double MaxFormFactor = 1.0/EnergyToWavelength(MaxE);
    MaxFormFactor += MaxFormFactor*0.02; // Add a 2% tolerance just to cover numerical errors.
    FormFactorData FormFactor( 0.0, MaxFormFactor, 1, 10000, FormFacFilename.c_str());*/
    
    
    std::string MuDataFilename;
    StringFromMap("AbsorptionData", InputData, MuDataFilename);
    
    double MinWavelength = EnergyToWavelength(MaxE); MinWavelength -= MinWavelength*0.02;
    double MaxWavelength = EnergyToWavelength(MinE); MaxWavelength += MaxWavelength*0.02;
    
    AbsorbCoeffData MuData( MinWavelength, MaxWavelength, 1, 5000, MuDataFilename.c_str());
  
    
    double a0 = 3.31; //lattice constant 2.87 angstroms (iron) 3.31 (Ta)
    DoubleFromMap("LatticeConstant", InputData, a0);

    Vector a1(1,0,0);
    Vector a2(0,1,0);
    Vector a3(0,0,1);

    a1 = a1*a0;
    a2 = a2*a0;
    a3 = a3*a0;

	
	std::string UnitCellFileName;
	StringFromMap("UnitCell", InputData, UnitCellFileName);
	LatticePlane::SetupUnitCell(UnitCellFileName.c_str(), a1, a2, a3, a0, MinE, MaxE, &MuData);
	

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

    //std::vector< hklSet > hklPlanes;
	std::vector< LatticePlane > hklPlanes;

    LatticePlane S(0,0,0,0);

    for( int n1 = 0; n1 <= indexMax; n1++)
    {
        if(n1 != 0)
        {
            S = LatticePlane(0 ,0 ,n1, 6);
            if( !isPlaneInList(S, &hklPlanes) )
            {
                hklPlanes.push_back(S);
                ScatterProbFile << 0  << 0  << n1 << "\t";
                BraggAngleFile << 0  << 0  << n1 << "\t";
                RockingCurveFile << 0  << 0  << n1 << "\t";
            }
            S = LatticePlane(0 ,n1,n1,12);
            if( !isPlaneInList(S, &hklPlanes) )
            {
                hklPlanes.push_back(S);
                ScatterProbFile << 0  << n1 << n1 << "\t";
                BraggAngleFile << 0  << n1 << n1 << "\t";
                RockingCurveFile << 0  << n1 << n1 << "\t";
            }
            S = LatticePlane(n1,n1,n1, 8);
            if( !isPlaneInList(S, &hklPlanes) )
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
                S = LatticePlane(0 ,n1,n2,24);
                if( !isPlaneInList(S, &hklPlanes) )
                {
                    hklPlanes.push_back(S);
                    ScatterProbFile << 0  << n1 << n2 << "\t";
                    BraggAngleFile << 0  << n1 << n2 << "\t";
                    RockingCurveFile << 0  << n1 << n2 << "\t";
                }
                S = LatticePlane(n1,n1,n2,24);
                if( !isPlaneInList(S, &hklPlanes) )
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
                    S = LatticePlane(n1,n2,n3,48);
                    if( !isPlaneInList(S, &hklPlanes) )
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

    ScatterProbFile << endl;
    BraggAngleFile << endl;
    RockingCurveFile << endl;
    
    double Temperature  = 300.0; //300K
    double DebyeTemperature = 240.0; //240K for Ta
    double mass_amu = 180.948; //180.948 amu for Ta
	
	DoubleFromMap("Temperature", InputData, Temperature);
	DoubleFromMap("DebyeTemperature", InputData, DebyeTemperature);
	DoubleFromMap("AtomicMass", InputData, mass_amu);
	
    
    double DebyeWallerPreFactor = LatticePlane::CalculateDebyeWallerPreFactor( Temperature, DebyeTemperature, mass_amu);
    cout << "DebyeWallerPreFactor:\t" << DebyeWallerPreFactor << endl;
	LatticePlane::SetDebyeWallerPreFactor( DebyeWallerPreFactor );
	
		
    for( int EnergyTick = 0; EnergyTick <= nEPoints; EnergyTick++)
    {
        double Energy = MinE + DeltaE*EnergyTick;

        ScatterProbFile << Energy << "\t";
        BraggAngleFile << Energy << "\t";
        RockingCurveFile << Energy << "\t";

        //float Sum = 0.0f;

        for( int i=0; i<int(hklPlanes.size()); i++)
        {
            //LatticePlane Plane( hklPlanes[i].h,  hklPlanes[i].k,  hklPlanes[i].l, hklPlanes[i].M);
			LatticePlane Plane = hklPlanes[i];
            BraggAngleFile << Plane.FindBraggReflectionAngle( EnergyToWavelength(Energy) ) << "\t";
            float I = Plane.CalculatePowderScatter( EnergyToWavelength(Energy) );
            //Sum += I;
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
    
    
    double SourceDivergence = 0.0;
    DoubleFromMap("SourceDivergence", InputData, SourceDivergence);
    
    
    float MinTheta = 0.0;
    float MaxTheta  = Deg2Rad( +SourceDivergence );
    float MinPhi = 0.0;
    float MaxPhi  = 2.0*PI;
    
    float MinCosTheta = cos( MinTheta );
    float MaxCosTheta = cos( MaxTheta );
    
    if (MinCosTheta > MaxCosTheta)
    {
        swap( MinCosTheta, MaxCosTheta);
    }
    
    if (MinPhi > MaxPhi)
    {
        swap( MinPhi, MaxPhi);
    }
    
    
    int NumThetaSteps = 1;
    int NumPhiSteps = 1;
    IntFromMap( "NumThetaSteps", InputData, NumThetaSteps);
    IntFromMap( "NumPhiSteps", InputData, NumPhiSteps);
    
    double deltaCosTheta = (MaxCosTheta-MinCosTheta)/double(NumThetaSteps);
    double deltaPhi = (MaxPhi-MinPhi)/double(NumPhiSteps);
  
    
    Vector Source( -3.0f, 0.0f, 5.0);    
    VectorFromMap("Source", InputData, Source);
    
    double MinBraggAngle = 10;
    double MaxBraggAngle = -1;

    
    Vector SourceToOrigin = -1.0*Source;
        
    for( int CosThetaStep = 0; CosThetaStep <= NumThetaSteps; CosThetaStep ++)
    {
        double CosTheta = MinCosTheta + double(CosThetaStep)*deltaCosTheta;
                
        for( int PhiStep = 0; PhiStep <= NumPhiSteps; PhiStep ++)
        {
            double Phi = MinPhi + double(PhiStep)*deltaPhi;
                                    
            Vector SourceToCrystal( acos(CosTheta), Phi, true, true);
            SourceToCrystal = TransformToNewFrame(SourceToCrystal, SourceToOrigin.GetTheta(), SourceToOrigin.GetPhi());
            SourceToCrystal = SourceToCrystal.Normalized();            
            
            Vector CrystalIntersection(Source.x - (Source.z/SourceToCrystal.z)*SourceToCrystal.x,
                                       Source.y - (Source.z/SourceToCrystal.z)*SourceToCrystal.y,
                                       0.0f); //intersection of ray with z=0 plane
            
            for(int iCCDCorner = 0; iCCDCorner < 4; iCCDCorner ++)
            {
                Vector CrystalToCCD = CCDCorners[iCCDCorner] - CrystalIntersection;
                CrystalToCCD = CrystalToCCD.Normalized();
                double BraggAngle = 0.5*acos( CrystalToCCD.Dot(SourceToCrystal) );//0.5 as the calculation gives scattering angle;
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
    
    cout << "MinAngle:\t" << Rad2Deg(MinBraggAngle) << "\tMaxAngle:\t" << Rad2Deg(MaxBraggAngle) << "\t(Degrees)" << endl;
    
    ofstream BraggAngleLimitsFile;
    BraggAngleLimitsFile.open( "BraggAngleLimits.txt" );
    
    //keep the first part as one word!
    BraggAngleLimitsFile << "MinimumBraggAngle:\t" << MinBraggAngle << endl;
    BraggAngleLimitsFile << "MaximumBraggAngle:\t" << MaxBraggAngle << endl;
    
    BraggAngleLimitsFile.close();
    
    cout << endl;
    cout << "///////////////////////////////////////////////////////////" << endl;
    cout << "//Finished running Diffraction Probabiities Pre-Processor//" << endl;
    cout << "///////////////////////////////////////////////////////////" << endl;
    cout << endl;
    
    
    return 0;

}
