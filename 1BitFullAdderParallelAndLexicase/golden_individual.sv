module individual(input wire a, input wire b, input wire ci, output wire sum, output wire co);
assign sum = !(((ci^b)^!(a)));
assign co = ((((a&ci)|b))&(((ci|a)|(a&ci))));

// Karnaugh map solution:
//assign sum = (a^b^ci);
//assign co  = ((b&ci)|(a&ci)|(a&b));
endmodule
