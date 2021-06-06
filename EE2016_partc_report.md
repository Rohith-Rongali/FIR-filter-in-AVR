[TOC]

# 						EE2016 part-c assignment report							    

#### 1. FIR Filters

FIR-"Finite Impulse Response "  filters, one of the two primary filters are those which have a finite impulse response. 

The impulse response is characterised by a set of filter coeffecients/tap weights. I have implemented a 5-tap and 32-tap bandpass FIR filter which are characterised by 5 and 32  filter coefficients respectively.

For a given sequence of inputs "x[n]", the output "y[n]" from the filter(with impulse response "h[n]") is given by the convolution of "x[n]" and "h[n]".
$$
y[n]=x[n]*h[n]
$$

$$
y[n]=x[n].a(0)+x[n-1].a(1)+....+x[n-N+1].a(N-1)
$$

a(0),a(1).... are filter coefficients and N is the no. of taps.

The number of Taps usually indicates the memory required, the amount of calculations to done while implementing and most importantly it means we will get a better passband.

DTFT is applied to  the outputs(y[n]) obtained from the AVR program which will allow us to see the working of the filter.
$$
y[n]------> Y(jw)
$$




#### 2. FIXED POINT ARTHIMATIC

|  b7  |  b6  |  b5  |  b4   |   b3   |   b2   |   b1   |   b0   |
| :--: | :--: | :--: | :---: | :----: | :----: | :----: | :----: |
|  1   | 0.5  | 0.25 | 0.125 | 2^{-4} | 2^{-5} | 2^{-6} | 2^{-7} |

The above is the fixed point representation used by me for the 8-bit no.'s in AVR. This can represent  signed numbers in the range of [-1,1)  in two's complement form.

The filter coefficients and the input samples are represented by 8-bits in the above described notation so that we can use FMULS(fractional multiplication for signed numbers). This is ok because all the filter coefficients  and the inputs are in the range of [-1,1), except for the DC signal inputs which should be equal to one is represented by an approximate number= 0111 1111 in decimal is equal to 0.99 . 





#### 3. DTFT

```python
import numpy as np
from scipy.fft import fft
import matplotlib.pyplot as plt
def DTFT(output_list,title):
    y=np.array(output_list)
    N = len(y)  #no of otp samples
    T = 1.0 / 8000.0
    Y = fft(y)
    w = np.linspace(0.0, 1.0/(2.0*T), N//2)
    
    plt.plot(w, 2.0/N * np.abs(Y[0:N//2]))
    plt.xlabel("frequency")
    plt.ylabel("magnitude")
    plt.title(title)
    plt.grid()

```





#### 4. IMPLEMENTATION DETAILS

##### a.  Inputs and output storage locations,  coefficients and circular buffer locations

I've considered 1000 samples for each of the three cases which I've inputted in 4 batches of nearly 250 each.(given the limitations of the memory size)

For maintaining continuity between each of the batches , last four samples from previous batch are taken to really represent the convolution between the 1000 inputs and the coefficients.

| input batch |                        no. of inputs                         |
| :---------: | :----------------------------------------------------------: |
|   Batch1    |               4 zeroes+250samples = 254inputs                |
|   Batch2    |    4 last samples from batch1  + 250 samples = 254 inputs    |
|   Batch3    |    4 last samples from batch2  + 250 samples = 254 inputs    |
|   Batch4    | 4 last samples from batch3  + 250 samples +4 zeroes  = 254 inputs |



Shows the input in 32-tap filter program<img src="Pictures\Annotation 2020-11-09 150936.png">

The input is given as in the the above image, the first few numbers(5 or 32) are the coefficients and the rest are inputs.

###### For 5-tap filter : 

The image shows the addresses for 5-tap coefficients and input buffer.<img src="Pictures\Annotation 2020-11-05 20003242.png">



Location to store output of 5-tap<img src="Pictures\Annotation 2020-11-05 20003242dsf.png">





values highlighted in yellow : coefficients,  green: circular buffer, blue: the outputs<img src="Pictures\Annotation 2020-11-05 20003242dsfsf.png">

The above image shows the memory during implementation of some batch of input for a 5-tap filter.





###### For 32-tap filter : 



The image shows the addresses for 32-tap coefficients and input circular buffer.<img src="Pictures\Annotation 2020-11-05 20003242dsfsfdfs.png">





Location to store output of 32-tap:<img src="Pictures\Annotation 2020-11-05 20003242dsfsy.png">





The output location<img src="Pictures\Annotation 2020-11-09 153911sfd.png">





yellow: buffer location (when loaded with dc inputs),  green: coefficients storage space<img src="Pictures\Annotation 2020-11-09 153911sfdabcd.png">



In both the programs, I used *X-pointer* for retrieving coefficient values while I used *Y-pointer* for 2 purposes  i.e.  for accessing output locations as well as for working with the buffer. For doing so I used registers "R11,R12" for maintaining the address of *output locations*, and register "R19" for the *buffer location*(the buffer location was defined with addresses under 1 byte).



##### b.  Circular buffer

The circular buffer is an array of twice the size of the values required so it will be 10 bytes for 5-tap and 64 bytes for 32-tap, the first half and the second half are identical. So that we can do seamless  incrementing(5 or 32-times) starting from a location in first half  and read the entire buffer, similarly while decrementing from a location in second half. The  working is more clearly illustrated  below:

|        |        |        |        |        |        |        |        |        |        |
| :----- | :----- | :----: | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| **a0** | a1     |   a2   | a3     | a4     | **a0** | a1     | a2     | a3     | a4     |
| a5     | **a1** |   a2   | a3     | a4     | a5     | **a1** | a2     | a3     | a4     |
| a5     | a6     | **a2** | a3     | a4     | a5     | a6     | **a2** | a3     | a4     |
| a5     | a6     |   a7   | **a3** | a4     | a5     | a6     | a7     | **a3** | a4     |
| a5     | a6     |   a7   | a8     | **a4** | a5     | a6     | a7     | a8     | **a4** |
|        |        |        |        |        |        |        |        |        |        |
| **a5** | a6     |   a7   | a8     | a9     | **a5** | a6     | a7     | a8     | a9     |
| a10    | **a6** |   a7   | a8     | a9     | a10    | **a6** | a7     | a8     | a9     |
| a10    | a11    | **a7** | a8     | a9     | a10    | a11    | **a7** | a8     | a9     |

Here *the values in bold are the oldest values in the buffer*, so if we maintain a pointer to address of one of the two old values say the old value in the right half, it will be enough, since the left one's address is always 5 less than the right one. And in both my programs I used R19 register to keep track of the old value in the right half.

*To update the buffer* replace the old values with the new one and increment the R19 to point to the next oldest value, also if we reach the end of buffer , instead of incrementing, load the address of the starting of right half of the buffer(as in the above figure it is moved from a4 to a5 when we reach the end).

*To read* the contents of the buffer:

 *from the oldest value to the newest*, we start from the oldest value in left half(whose address can be obtained by subtracting 5 from R19)   and go on incrementing for 5 times to get all the 5 values in the buffer.

 *from the newest value to the oldest*, we pre-decrement starting from the oldest value in the right half(pointed by R19)  for 5 times to get all the values in the buffer.



**The above explanation can be extended to buffer of any size and hence for the 32-tap case as well**.



##### c. Details of the multiplication 

The coefficients and inputs as mentioned before are 8 bit signed no.'s in range [-1,1), their multiplication(which will  still be in the range [-1,1) ) by "FMULS" gives the output of 16-bit in the registers R0,R1 where the extra byte(lowerbyte:R0) increases the precision of the decimal , R1 has the same fixed point notation mentioned before.

|  b7  |  b6  |  b5  |  b4   |   b3   |   b2   |   b1   |   b0   |
| :--: | :--: | :--: | :---: | :----: | :----: | :----: | :----: |
|  1   | 0.5  | 0.25 | 0.125 | 2^{-4} | 2^{-5} | 2^{-6} | 2^{-7} |

this is for the upper byte R1.



|   b7   |   b6   |   b5    |   b4    |   b3    |   b2    |   b1    |   b0    |
| :----: | :----: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: |
| 2^{-8} | 2^{-9} | 2^{-10} | 2^{-11} | 2^{-12} | 2^{-13} | 2^{-14} | 2^{-15} |

This is for the lower byte R0.



The result is accumulated in R4,R5 and if any carry results while addition is stored in R6. 

For 5-tap five such multiplication results are added and accumulated in each MAC cycle, similarly 32 for 32-tap. While adding, it will be like adding unsigned numbers. So, a register R3 is maintained to keep track of the number of negative values resulted from multiplication, and subtract R3 from R6 to get the correct result.



##### d. A brief of the program logic

 If 'N' inputs are given 'N-4' MAC cycles are needed and each MAC cycle is iterated 5 or 32 times.

X-pointer initially points to start of the coefficients buffer, and after one MAC cycle it points to the end .So, in next MAC cycle we keep decrementing to reach the start of the buffer,and the same thing repeats. So, the pointer moves from a(0) to a(4)  and  a(4) to a(0) in alternate cycles. Accordingly to this the pointer to  input circular buffer must move from newest value to oldest value and oldest value to newest value in alternate cycles. To identify this we maintain a register initially equal to zero and one's complement this register after each MAC cycle. Based on the value in this register the flow of input access and the coefficient access is decided at the start of each cycle.



R22 is used for this purpose of the identifying the flow of memory access<img src='Pictures\jdkfj.png'>







Below are the instructions at the end of each MAC cycle initialising the y pointer for next MAC cycle for input access and also complements the bits in R22 to change the flow of memory access in the next MAC cycle.

<img src="Pictures\sdfgsafsa.png" alt="sdfgsafsa"  />



Finally at the end of each MAC cycle, the output from R6,R5 are stored to locations in the registers R11,R12. And also the addresses in them are incremented for the outputs of the next MAC cycle.

<img src="Pictures\dnfkks.png">

#### 5.RESULTS

##### a. Input and output plot



input : blue, output: yellow<img src="Pictures\sine5.png">





input : blue, output: red<img src="Pictures\sine32.png">





input : blue, output: red<img src="Pictures\wgn5.png">





input : blue, output: red<img src="Pictures\wgn32.png">







input : blue, output: red<img src="Pictures\dcout5.png">





input : blue, output: red<img src="Pictures\dcout32_c.png">









##### b. Frequency response plots



###### 1.DC signal:



<img src="Pictures\fdcinp.png">





<img src="Pictures\dc5outf.png">





<img src="Pictures\dc32outf.png">







###### 2.SINE @ 1800Hz:

<img src="Pictures\sininpf.png">





<img src="Pictures\sine5otf.png">





<img src="Pictures\sine32otf.png">



###### 3. White noise:



<img src="Pictures\wgninpf.png">



<img src="Pictures\wgn5otf.png">

<img src="Pictures\wgn32otf.png">