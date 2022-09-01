// Definitions
`define TEST_COUNT 8
`define PERIOD 10

// Use parameter to pass population size to testbench
module testbench();
  // Use integer array to store fitness for all individuals
  integer fitness;

  // Define inputs and outputs to connect to the individuals
  reg clk;
  reg rst;
	reg a;
	reg b;
	reg ci;
	reg sum_expected;
	wire sum_current;
	reg co_expected;
	wire co_current;

	// Array of testvectors to store expected values
	reg [4:0] testvectors[0:(`TEST_COUNT-1)];
	integer vectornum = 0;
	
	logic [0:(`TEST_COUNT-1)] lexi_sum;
	logic [0:(`TEST_COUNT-1)] lexi_co;
	
	reg hit_sum;
	reg hit_co; 

	// Instantiate all the individuals
       individual dut_1(.a(a), .b(b), .ci(ci), .sum(sum_current), .co(co_current));
 

	// Create task to evaluate each testcase
	// Note that this evaluates the entire population for this testcase
	task testcase;
		input a_value, b_value, ci_value, sum_value, co_value;

		// Initialise inputs
		a=a_value;
		b=b_value;
		ci=ci_value;
		
		hit_sum = 0;
		hit_co = 0;

		#(`PERIOD/2) // Wait till negative clock edge to check signals
		
      		if(sum_current == sum_value) 
      		begin
        		fitness = fitness + 1;
        		hit_sum = 1;
		end
		
		if(co_current == co_value) 
      		begin
        		fitness = fitness + 1;
        		hit_co = 1;
		end
	endtask

	// We want to change our inputs on negedge so that we can run testcases at posedge
	always @(posedge clk) begin
		if (rst == 0) begin
			{a, b, ci, sum_expected, co_expected} = testvectors[vectornum];
			
			testcase(a, b, ci, sum_expected, co_expected);
			//$display("Vectornum: %d Hit_sum: %d, Hit_co: %d",vectornum, hit_sum, hit_co);
			
			
			lexi_sum[vectornum] = hit_sum;
			lexi_co[vectornum] = hit_co;
			
			vectornum = vectornum + 1;
	
		end

	end

  // We use a clock to load in our testvectors
  always begin
    #(`PERIOD/2) clk = ~clk;
  end

  // Run simulation
  initial begin
    // Set the clock high so we get posedges at 10,20,30 etc.
    //clk = 1;
    clk = 1; rst = 1;
    #(`PERIOD) rst = 0;

    // DEBUGGING - Dump to vcd file
    //$dumpfile("testbench_values.vcd");
    //$dumpvars(0,fulladder_tb);

    // Read in test vectors
    $readmemb("test_vectors.tv", testvectors);

    //$display("%0b",testvectors[4]);
    
    // Set all scores to zero
    fitness = 0;
   
     
    // Wait 8 clock cycles for test to complete
    #(`PERIOD*`TEST_COUNT*2);

    // Print all fitness scores to the console
    //$display("%0d",fitness);
    // Print all fitness scores to the console
    //$display("Lexi_sum: %b",lexi_sum);
    //$display("Lexi_co: %b",lexi_co);
    $display("%0b%0b",lexi_sum,lexi_co);

    $finish;
  end
 endmodule
 


