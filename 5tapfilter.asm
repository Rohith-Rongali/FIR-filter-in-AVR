;
; filter5.asm
;
; Created: 04-11-2020 11:09:49
; Author : HP
;

;LET'S SUPPOSE THAT THE ACCUMULATOR RESTS AT R4,R5. R3 CONTAINS THE '-2' COUNTER.

; Replace with your application code

; Replace with your application code






LDI ZL,LOW(NUM<<1);
LDI ZH,HIGH(NUM<<1);
LDI R24,0X05;  load the number of coefiicients

LDI XL,0X60;
LDI XH,0X00;location for coefficients 60-64

LDI YL,0XA0;
LDI YH,0X00; location to store outputs
MOV R11,YL;
MOV R12,YH;

LDI R19,0X76; ADDRESS FOR THE BUFFER


;Loading of coefficients into **X:0X60-0X64**
;R23 JUST USED AS A SPARE REGISTER, R24 COUNTER FOR COEFFICIENTS, FINALLY BECOMES ZERO.
CONTLD:LPM R23,Z+;
ST X+,R23;
DEC R24;
BREQ CONTLDOUT;
RJMP CONTLD;
CONTLDOUT:

MOV YL,R19;
LDI R24,5; AGAIN AS A COUNTER.
;LOADING VALUES INTO BUFFER "Y:0X76-0X7F".
;R23 JUST USED AS A SPARE.
CONTPUSH:LPM R23,Z+;
STD Y+5,R23;
ST Y+,R23;
DEC R24;
BREQ STOPPUSH;
RJMP CONTPUSH;
STOPPUSH:
;Y NOW POINTS TO OLD
MOV R19,YL;
LDI YH,0X00;


LDI R23,0XFE; NO OF OUTPUTS COUNTER

LDI R22,0X05; JUST SPARE
MOV R9,R22; $$$$ R9 I ALWAYS EQUAL TO NO. OF TAP, IT'S A CONSTANT.

LDI R22,0X00; IMPORTANT KEEP IT HERE ONLY(RECOGNISES DEC OR INC CYCLE)



SUB YL,R9;

;;;;;;STARTING OF THE BIGLOOP;;;;;;
 
CONTBIGLOOP:

LDI R21,0X00; R21 JUST A SPARE NO SIGNIFICANCE
MOV R3,R21; '-2' COUNTER  ;CLEARING ALL THE REGISTERS BEFORE MAC CYCLE
MOV R4,R21; LOWER BYTE
MOV R5,R21; UPPER BYTE
MOV R6,R21; CARRY HOLDER



LDI R24,0X05; NO OF MULTIPLICATIONS COUNTER(=NO OF COEFFECIENTS)


;;;;;;;;;;; MULTIPLY AND ACC;;;;;;;;;;MAC CYCLE


CONTMAIN:

SBRC R22,0;
RJMP INCDEC;
LD R16,-X; 
LD R17,Y+;
;ST Y+,R17;

RJMP DECINC;

INCDEC:
LD R16,X+;
LD R17,-Y;

DECINC:

FMULS R16,R17;

;ADDING R0,R1 TO R4,R5,R6 :--
;ADD R4,R0;
;BRCC NOCARRYINADDL;
;INC R5;
;BRCC IFCARRY;
;INC R6;
;IFCARRY:
;NOCARRYINADDL:ADD R5,R1;
;BRCC NOCARRYINADDU;
;INC R6;
;NOCARRYINADDU:

;ADDING R0,R1 TO R4,R5,R6 :--
ADD R4,R0;
ADC R5,R1;
BRCC NOCARRY;
INC R6;
NOCARRY:




LSL R1; UPPER BIT OF THE MULTIPLICATION 
BRCC NOTNEG1;
INC R3;
NOTNEG1:
DEC R24;
BREQ STOP;

RJMP CONTMAIN;

;;;;;;;;;;; MULTIPLY AND ACC;;;;;;;;;;;;

STOP:



;STORE THE OUTPUTS AT THE DESIRED LOCATIONS--
SUB R6,R3;
MOV YL,R11;
MOV YH,R12;
ST Y+,R6;
ST Y+,R5;
MOV R11,YL;
MOV R12,YH;

;REPLACING THE OLD WITH NEW
MOV YL,R19;
LDI YH,0X00;
LPM R17,Z+;
ST Y,R17; REPLACING THE OLD ONE WITH THE NEW
SUB YL,R9; R9=5 YL IS THE LEFT OLD.
ST Y,R17;


;INC TO GET THE NEXT OLD
INC R19;
SBRC R19,7;
LDI R19,0X7B;



SBRC R22,0;
RJMP PDECINC;
MOV YL,R19;   MAKING Y POINT TO THE RIGHT OLD.
RJMP PINCDEC;

PDECINC:
INC YL; MAKING Y POINT TO THE LEFT OLD/Right OLD making it feasible to increment in the next cycle
PINCDEC:



COM R22; complementing R22 useful for changing the memory access scheme in the next cycle.

DEC R23;
BREQ FULLSTOP;

RJMP CONTBIGLOOP;

FULLSTOP:


NOP;

NUM: .db 0XCE,0X0C,0X58,0X0C,0XCE,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x7F,0x0,0x0,0x0,0x0; COEFFICIENTS AND INPUTS