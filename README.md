# Composition-of-Active-and-Passive-Circuits

### Description
This collection of scripts was created for the scope of an assignment in a Composition-of-Active-and-Passive-Circuits cource of the Aristotle University of Thessaloniki during the 2018-19 academic year. The main purpose of these projects is to demonstrate the implementation of a series of different active circuits - analog filters, satisifying a set of specifications. For this purpose important intel such as FFT Analysis, Transient Analysis, Bode Diagrams etc. is provided.

### Implementations - Specifications

#### Low Pass Filter - Inverse Chebyshev
 - f_p  = 5.5 KHz
 - f_s = 10.45 KHz
 - a_max = 0.3 dB
 - a_min = 25 dB
 - LowFrequencyGain = 5 dB
 
![alt-text-1](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/LP/plots/T_LP.png "Low Pass Filter")
![alt-text-2](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/LP/plots/finalCircuit.jpg "Low Pass Circuit")


#### Band Pass Filter - Butterworth
 - f_0 = 1.2 KHz
 - f_1 = 675 Hz
 - f_2 = 2.1333 KHz
 - f_3 = 326.54 Hz
 - f_4 = 4409.87 Hz
 - a_max = 0.5278 dB
 - a_min = 26.7778 dB
 - LowFrequencyGain = 10 dB
 
![alt-text-3](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/BP/plots/T_BP.png "Band Pass Filter")
![alt-text-4](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/BP/plots/FinalCircuit.jpg "Band Pass Circuit")


#### Band Elimination Filter - Inverse Chebyshev
 - f_0 = 2.4 KHz
 - f_1 = 1.75 KHz
 - f_2 = 3291.4286 Hz
 - f_3 = 2075.1088 Hz
 - f_4 = 2775.7581 Hz
 - a_max = 0.5556 dB
 - a_min = 30.7778 dB
 - LowFrequencyGain = 10 dB
 
![alt-text-5](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/BE/plots/T_BE.png "Band Elimination Filter")
![alt-text-6](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/BE/plots/finalCircuit.jpg "Band Elimination Circuit")


#### High Pass Filter - Inverse Chebyshev
 - f_p  = 3 KHz
 - f_s = 1.667 KHz
 - a_max = 0.5278 dB
 - a_min = 27.2222 dB
 - LowFrequencyGain = 0 dB
 
![alt-text-7](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/HP/plots/T_HP.png "High Pass Filter")
![alt-text-8](https://github.com/kosletr/Composition-of-Active-and-Passive-Circuits/blob/master/HP/plots/finalCircuit.jpg "High Pass Circuit")

### Setup
The provided code was created using MATLAB R2019a and MULTISIM 14.0, however older software versions should work fine. All of the .m scripts provided, are  commented for higher readability and maintenance.
