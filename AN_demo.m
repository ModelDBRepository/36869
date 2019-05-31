%  AN_Demo.m
%
%
%	 							Dual-AGC Model for Auditory Nerve Adaptation
%
%													D.C. Mountain
%												Boston University
%											Hearing Research Center
%
%	            ____________     ___________     __________     __________
%	           |            |   |           |   |          |   |           |
%	input ---->| Ca channel |-->|  Vesicle  |-->| Synaptic |-->| Receptor  |--> output
%	           |            |   | Depletion |   |   Cleft  |   | Depletion |
%	            ------------     -----------     ----------     ----------
%
%													MODEL STRUCTURE
% REFERENCES
%
% Deligeorges, S. and Mountain D.C. (1997) A Model for Periodicity Coding in the Auditory System.
%     In: Computational Neuroscience: Trends in Research, 1997, J.M Bower ed., Plenum, New York pp 609-616. 
% Singh, S. and Mountain, D.C. (1997) A Model for Duration Coding in the Inferior Colliculus. 
%     In: Computational Neuroscience: Trends in Research, 1997, J.M Bower ed., Plenum, New York pp 497-503. 



% These parameters are to simulate a high spontaneous rate fiber

Ks = 0.8;		%(1/mV)  the voltage sensitivity of the Ca channel
V0 = 2.1;		%(mV)  voltage offset of the Ca channel re IHC resting potential
        			%these parameters give the synapse a dynamic range of ~5 mV
        
K3 = 375.0;     %{vesicles/s}
K4 = 125.0;     %{vesicles/s}
        
TauCleft =0.0001;  %(s) synaptic cleft time constant
        
Kr = 33.0;		%{1/s}
Kp = 8.16;		%{receptors/s}
N0 = 36.09;		%{receptors}
        
% Prompt the user for the input level        
Level = input('Input Level from IHC (mV)  ==>');
        
deltaT = 1/50000;		% Time step size in seconds
FrameSize = 5000;		% This combination creates a 100 ms frame

% The following constants are for the digital filters used in the model.
% They were derived from the continuous-time equations using the bilinear transform.

Ka1    = 1/(2.0/(K4*deltaT) + 1);
Kb1    = 2.0/(K4*deltaT) - 1;

Ka2    = 1/(2.0*TauCleft/deltaT + 1);
Kb2    = 2.0*TauCleft/deltaT - 1;

Ka3    = 1/(2.0/(Kp*deltaT) + 1);
Kb3    = 2.0/(Kp*deltaT) - 1;

Time = 1000*deltaT*[0:FrameSize-1]; % Time array for plotting in ms

Nstart = 500;								% input turns on after 500 time steps
Nstop = FrameSize-2000;					% input turns off 2000 time steps from end of frame

DataIn 	= zeros(1,FrameSize);								% contains IHC receptor potential waveform
DataIn(Nstart:Nstop)	= Level*ones(1, Nstop-Nstart+1);	% load the array with a pulse to represent a tone burst

DataOut 	= zeros(1,FrameSize);		% will contain output data (spikes/s)

%	Start the time loop, but first set the initial conditions

B     = 0.03;
Filt1 = 0.0925;          %{this is equal to 1-A-B}
Filt2 = 0.03;
R     = 32.44;
Filt3 = N0 - 32.1144;      %{this is equal to N0-N}

for i = 1 : FrameSize					%begin integration loop

               Vm  = DataIn(i);		%IHC membrane potential (Vrest=0)
               if  Vm > 18.0 
                  Vm = 18.0; 			%to avoid F.P. overflow
                  end;

% Ca Channel open probability
               P        = 1/(1+exp(2*Ks*(-Vm+V0)));
               
% Vesicle depletion
               Bold     = B;
               B        = (1 - Filt1)*P;	%vesicle release rate
               OldFilt1 = Filt1;
               Filt1    = Ka1*((K3/K4)*(B+Bold)+Kb1*OldFilt1);	% low-pass filter

% Synaptic Cleft
               OldFilt2 = Filt2;
               Filt2    = Ka2*(B+Bold+Kb2*OldFilt2);	%neurotransmitter concentration in cleft

% Receptor depletion
               Rold     = R;
               R        = Filt2*Kr*(N0-Filt3);	%firing rate is proportional to number of active receptors binding neurotransmitter
               OldFilt3 = Filt3;
               Filt3    = Ka3*((1/Kp)*(R+Rold)+Kb3*OldFilt3);	% low-pass filter

% Load the output file
               DataOut(i) = R;


end;  %integration loop

%Generate a graph of the instantaneous firing rate

plot(Time ,DataOut);
axis([0 max(Time) 0 1000]);
xlabel('Time (ms)','FontSize',12,'FontWeight','bold');
ylabel('Firing Rate (spikes/s)','FontSize',12,'FontWeight','bold');
title('Auditory-Nerve Adaptation','FontSize',12,'FontWeight','bold');
shg;