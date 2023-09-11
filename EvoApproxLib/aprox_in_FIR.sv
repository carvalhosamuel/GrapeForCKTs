`timescale 1ns / 1ps

// 16-Bit Approximate adders 
`include "16bit_adders/add16se_2T2.v"  
`include "16bit_adders/add16se_32T.v"
`include "16bit_adders/add16se_2TN.v"  
`include "16bit_adders/add16se_334.v"
`include "16bit_adders/add16se_2U6.v"  
`include "16bit_adders/add16se_349.v"
`include "16bit_adders/add16se_2UB.v"  
`include "16bit_adders/add16se_35M.v"
`include "16bit_adders/add16se_2UY.v"  
`include "16bit_adders/add16se_36D.v"
`include "16bit_adders/add16se_2X2.v"  
`include "16bit_adders/add16se_3BA.v"
`include "16bit_adders/add16se_2YM.v"  
`include "16bit_adders/add16se_3BD.v"
`include "16bit_adders/add16se_32R.v"

// 8-Bit Approximate multipliers 
`include "8bit_multipliers/mul8s_1L12.v"  
`include "8bit_multipliers/mul8s_1KV6.v"  
`include "8bit_multipliers/mul8s_1KVA.v"  
`include "8bit_multipliers/mul8s_1KVQ.v"  
`include "8bit_multipliers/mul8s_1L2N.v"
`include "8bit_multipliers/mul8s_1KV8.v"  
`include "8bit_multipliers/mul8s_1KVM.v"  
`include "8bit_multipliers/mul8s_1KX5.v"  
`include "8bit_multipliers/mul8s_1L2J.v"
`include "8bit_multipliers/mul8s_1KV9.v"  
`include "8bit_multipliers/mul8s_1KVP.v"  
`include "8bit_multipliers/mul8s_1KXF.v"  
`include "8bit_multipliers/mul8s_1L2L.v"


// This is the testbench for the "golden individual", with the coefficients already given
// It implements a moving average filter with all four coefficients set to 8'h20, as seen below
// This will be the target for GE to achieve

// Coefficients that will be found by the grammar

// Flip flop definition

module DFF(clk, reset, data_in, data_delayed);
parameter N = 8;
input clk, reset;
input signed [N-1:0] data_in;
output reg signed [N-1:0] data_delayed; 

always@(posedge clk, posedge reset)
begin
    if (reset)
    data_delayed <= 0;
    else
    data_delayed <= data_in; 
    
end

endmodule


// Testbench definition

module FIR_TB;

parameter coef_0 =  8'h20;
parameter coef_1 =  8'h20;
parameter coef_2 =  8'h20;
parameter coef_3 =  8'h20; 
 
parameter N = 8;

reg clk, reset;
reg [N-1:0] data_in;
wire [15:0] data_out;
reg [N-1:0] RAMM [31:0];  
reg [4:0] Address; 
reg [7:0] b0 =  coef_0; 
reg [7:0] b1 =  coef_1; 
reg [7:0] b2 =  coef_2; 
reg [7:0] b3 =  coef_3;

FIR_Filter inst0(clk, reset, b0, b1, b2, b3, data_in, data_out);

// input sine wave data
initial
// Samuel: Added 0,31 at the end to supress warning of IEEE standard
$readmemb("signal.data", RAMM,0,31);

// create a clock
initial 
clk = 0;
always 
#10 clk = ~ clk;  

// Read RAMM data and give to design
always@(negedge clk)
    data_in <= $signed(RAMM[Address]);  	
// Address counter
initial
Address = 0; 
always@(posedge clk)
begin
    if (Address == 31)
      begin
        Address = 0;
        $finish;
      end
    else
        Address = Address + 1; 
end     

endmodule
 
 
// Individual definition 
 
module FIR_Filter(clk, reset, b0, b1, b2, b3, data_in, data_out);

parameter N = 8;

input clk, reset;
input signed [7:0] b0, b1, b2, b3;
input signed [N-1:0] data_in;
output reg signed [15:0] data_out; 

// coefficients defination
// Moving Average Filter, 3rd order
// four coefficients; 1/(order+1) = 1/4 = 0.25 
// 0.25 x 128(scaling factor) = 32 = 6'b100000

wire signed [N-1:0] x1, x2, x3; 

// Create delays i.e x[n-1], x[n-2], .. x[n-N]
// Instantiate D Flip Flops
// Samuel: changed 0 to 1'b0 on 
DFF DFF0(clk, 1'b0, data_in, x1); // x[n-1]
DFF DFF1(clk, 1'b0, x1, x2);      // x[x[n-2]]
DFF DFF2(clk, 1'b0, x2, x3);      // x[n-3] 
//DFF DFF3(clk, 1'b0, x3, x4);      // x[n-4]

//  Multiplication
wire signed [15:0] Mul0, Mul1, Mul2, Mul3,Mul4;  


// Instantiating approximate multipliers
mul8s_1KV6 mul0 (data_in,b0,Mul0);
mul8s_1KV6 mul1 (x1,b1,Mul1);
mul8s_1KV6 mul2 (x2,b2,Mul2);
mul8s_1KV6 mul3 (x3,b3,Mul3);


//assign Mul0 = data_in * b0; 
//assign Mul1 = x1 * b1;  
//assign Mul2 = x2 * b2;  
//assign Mul3 = x3 * b3;  
//assign Mul4 = x4 * b4;
 
// Addition operation
wire signed [15:0] Add_final; 
wire signed [15:0] add1;
wire signed [15:0] add2;
wire signed [15:0] add3;


// Instantiate the approximate design and connect to testbench variables  
add16se_2TN  fa1 ( Mul0, Mul1, add1);  
add16se_2TN  fa2 ( add1, Mul2, add2);
add16se_2TN  fa3 ( add2, Mul3, add3);  

assign Add_final = Mul0 + Mul1 + Mul2 + Mul3; 
//assign Add_final = add3;

// Final calculation to output 
always@(posedge clk)
// Samuel: Added begin/end as multiple statements needed with $monitor included to get output
begin
	data_out <= Add_final; 
	$monitor("%b", data_out);
	//$display("%b", data_in);
end
endmodule
