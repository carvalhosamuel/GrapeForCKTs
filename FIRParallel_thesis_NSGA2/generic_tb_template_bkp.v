// This is the testbench for the "golden individual", with the coefficients already given
// It implements a moving average filter with all four coefficients set to 8'h20, as seen below
// This will be the target for GE to achieve

// Coefficients that will be found by the grammar

// PARAMETERS:
// N = Number of bits
// S = Number of samples in signal


`timescale 1ns / 1ps


// Flip flop definition


module DFF(clk, data_in, data_delayed);
parameter N = 32;
input clk;
input [N-1:0] data_in;
output reg [N-1:0] data_delayed=0; 

always@(posedge clk)
begin
    /*if (reset)
    data_delayed <= 0;
    else*/
    data_delayed <= data_in; 
    
end

endmodule


// Testbench definition

module FIR_TB;

//COEFFS HERE

//DYNAMIC1

//END DYNAMIC1

parameter N = 32;
parameter S = 128;

reg clk, en;
reg [N-1:0] data_in;
wire [N-1:0] data_out;
reg [N-1:0] RAMM [S-1:0];  
reg [7:0] Address;


//DYNAMIC2

//END DYNAMIC2

//DYNAMIC3

//END DYNAMIC3



// input sine wave data
initial
// Samuel: Added 0,31 at the end to supress warning of IEEE standard
$readmemb("signal_noisy_128_sine.data", RAMM,0,S-1);

// create a clock
initial 
clk = 0;
always 
#10 clk = ~ clk;  

// Read RAMM data and give to design
always@(posedge clk)
begin
    data_in <= RAMM[Address];  
    //Ideally enable should only be on after "number of taps" clock cycles, but we will workaround it on the python fitness function
    en <= 0;
    $monitor("%b", data_out);
end
// Address counter

initial
Address = 0; 
always@(posedge clk)
begin
    if (Address == S-1)
      begin
        Address = 0;
        $finish;
      end
    else
        Address = Address + 1; 
end     

endmodule
 
 
// Individual definition 
 
//DYNAMIC4

//END DYNAMIC4

parameter N = 32;

//DYNAMIC5

//END DYNAMIC5

//DYNAMIC6

//END DYNAMIC6

input clk, en;
input [N-1:0] data_in;
output reg [N-1:0] data_out; 



// Create delays i.e x[n-1], x[n-2], .. x[n-N]
// Instantiate D Flip Flops

//DYNAMIC7

//END DYNAMIC7

//Multiplication

//DYNAMIC8 

//END DYNAMIC 8   
 
// Addition operation
wire [N-1:0] Add_final;
 
//DYNAMIC9 

//END DYNAMIC 9 

// Final calculation to output 
always@(posedge clk)
begin
if(en)
	data_out <='hz;
else if(!en)
begin
	data_out <= Add_final;
	//$monitor("%b", data_out);
end	 
else
begin
	data_out <=data_out;
	//$monitor("%b", data_out);
end
end
endmodule
