'''
Processes the data from DiffractionProb.


'''
import sys


print "//////////////////////////////"
print "//Processing Scattering Data//"
print "//////////////////////////////"


bForceDiffraction = False;

if(len(sys.argv) == 2):
	if( sys.argv[1] == "-ForceDiffraction" ):
		print "Diffraction forcing is on!"
		bForceDiffraction = True;



#0.610865238; #35 degrees
#0.785398163; #45 degrees
#0.34906585; #20 degrees
#0.872664626; #50 degrees
MinAngle = 0.717;
MaxAngle = 0.73;

AngleLimitData = [line.split('\t') for line in open('BraggAngleLimits.txt')];

MinAngle = float(AngleLimitData[0][1]);
MaxAngle = float(AngleLimitData[1][1]);

print "Angle Ranges:\t", MinAngle, "-->\t", MaxAngle

AngleData = [line.split('\t') for line in open('BraggAngle.txt')];
for l in AngleData:
	del l[-1]; #need to get rid of newline character at the end
	
ProbData = [line.split('\t') for line in open('ScatterProb.txt')];
for l in ProbData:
	del l[-1];

RockingCurveData = [line.split('\t') for line in open('RockingCurve.txt')];
for l in RockingCurveData:
	del l[-1];

KeepPlane = [];

for i in range( 1, len(AngleData[0]) ): 
	KeepPlane.append(False);


for i in range( 1, len(AngleData[0]) ): #cycle from 1 to last part of a line
	#loop down each angledata line
	for AngleLine, ProbLine in zip(AngleData[1:],ProbData[1:]):
		if MaxAngle > float(AngleLine[i]) > MinAngle: # and float(ProbLine[i]) > 1E-5:			
			KeepPlane[i-1] = True;

RemainingPlanes = [];


for bShouldKeep,PlaneName in zip(KeepPlane,AngleData[0][1:] ):
	if(bShouldKeep):
		RemainingPlanes.append( PlaneName );
		
print "Remaining Planes:", RemainingPlanes;
print KeepPlane;


Output = open("ProcessedScatteringData.txt" ,'w');

Output.write( "Energy\t" );

for PlaneName in RemainingPlanes:
	Output.write( PlaneName + "\t\t\t");

Output.write("Sum\n");

OutputArray = []


#Go down by rows (scan in energy) followed by columns (lattice planes)
for angleline, probline, rockingcurveline in zip(AngleData[1:],ProbData[1:],RockingCurveData[1:]):
	OutputLine = [] 
	OutputLine.append(angleline[0])
	#Output.write(angleline[0] + "\t"); #print energy for this row
	ProbSum = 0.0;
	for prob,angle in zip(probline[1:],angleline[1:]):
		#if MaxAngle > float(angle) > MinAngle: #Only include reflections that we care about
		ProbSum += float(prob); #sum up all probabilities across row
			
	for angle,prob,rockingcurve,i in zip(angleline[1:],probline[1:],rockingcurveline[1:],range(len(probline[1:]))):
		
		if( KeepPlane[i] ): #only print this column if we're keeping this plane			
			if( ProbSum != 0.0 and MaxAngle > float(angle) > MinAngle ):   
				OutputLine.append(angle)
				OutputLine.append(float(prob)/ProbSum)   
				OutputLine.append(rockingcurve)		  
				#Output.write( angle + "\t" + str(float(prob)/ProbSum) + "\t" + rockingcurve + "\t");
			else:
				#make probability zero if we don't care about the angle
				OutputLine.append(angle)
				OutputLine.append(0)
				OutputLine.append(0)
				#Output.write( angle + "\t" + str(0)				   + "\t" + str(0)	   + "\t");
				
	OutputLine.append(ProbSum)
	OutputArray.append(OutputLine)
	#Output.write( str(ProbSum) + "\n");
	

if bForceDiffraction:
	MaxProb = -1.0
	
	for line in OutputArray:
		if float(line[-1]) > MaxProb :
			MaxProb = line[-1]
	
	for line in OutputArray:
		line[-1] = line[-1]/MaxProb		
	

for line in OutputArray:
	for element in line[:-1]:
		Output.write(str(element) + "\t")
		
	Output.write(str(line[-1]) + "\n")
		
Output.close();

Output = open("DiffractionPeakLimits.txt" ,'w');

#code to determine diffraction peaks/energy ranges
#scan along columns. If the plane is being kept then scan in energy to find the range of diffracted energies that are expected to hit the CCD
for bShouldKeep,PlaneName,index in zip(KeepPlane,AngleData[0][1:],range(len(KeepPlane)) ): #scan across planes
	if(bShouldKeep):
		MinEnergy = 1000
		MaxEnergy = -1
		for angleline in AngleData[1:]: #scan in energy
			if (MaxAngle > float(angleline[index+1]) > MinAngle) :
				if float(angleline[0]) < MinEnergy:
					MinEnergy = float(angleline[0])
				if float(angleline[0]) > MaxEnergy:
					MaxEnergy = float(angleline[0])
		print PlaneName, ":\t", MinEnergy, "--->", MaxEnergy
		Output.write( PlaneName + "\t" + str(MinEnergy) + "\t" + str(MaxEnergy) + "\n" )


Output.close();

print "///////////////////////////////////////"
print "//Finished Processing Scattering Data//"
print "///////////////////////////////////////\n"










