#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <fstream>

using namespace std;

#define MemSize 1000 // memory size, in reality, the memory size should be 2^32, but for this lab, for the space resaon, we keep it as this large number, but the memory is still 32-bit addressable.
#define add 1
#define sub 2
#define and 3
#define or  4
#define xor 5


struct IFStruct
{
	bitset<32> PC;
	bool nop;
};

struct IDStruct
{
	bitset<32> Instr;
	bool nop;
};

struct EXStruct
{
	bitset<32> Read_data1;
	bitset<32> Read_data2;
	bitset<16> Imm;
	bitset<5> Rs;
	bitset<5> Rt;
	bitset<5> Wrt_reg_addr;
	bool is_I_type;
	bool rd_mem;
	bool wrt_mem;
	bool alu_op; // 1 for addu, lw, sw, 0 for subu
	bool wrt_enable;
	bool nop;
};

struct MEMStruct
{
	bitset<32> ALUresult;
	bitset<32> Store_data;
	bitset<5> Rs;
	bitset<5> Rt;
	bitset<5> Wrt_reg_addr;
	bool rd_mem;
	bool wrt_mem;
	bool wrt_enable;
	bool nop;
};

struct WBStruct
{
	bitset<32> Wrt_data;
	bitset<5> Rs;
	bitset<5> Rt;
	bitset<5> Wrt_reg_addr;
	bool wrt_enable;
	bool nop;
};

struct stateStruct
{
	IFStruct IF;
	IDStruct ID;
	EXStruct EX;
	MEMStruct MEM;
	WBStruct WB;
};

struct forwardingStruct {
    // 2-bit control signal for ALU input A:
    // 0: No forwarding, 1: forward from EX/MEM, 2: forward from MEM/WB.
    unsigned int forwardA = 0;

    // 2-bit control signal for ALU input B:
    // 0: No forwarding, 1: forward from EX/MEM, 2: forward from MEM/WB.
    unsigned int forwardB = 0;
};


class InsMem
{
public:
	string id, ioDir;
	int size;
	InsMem(string name, string ioDir)
	{
		id = name;
		IMem.resize(MemSize);
		ifstream imem;
		string line;
		size = 0;
		imem.open(ioDir + "\\imem.txt");
		if (imem.is_open())
		{
			while (getline(imem, line))
			{
				IMem[size] = bitset<8>(line);
				size++;
			}
		}
		else
			cout << "Unable to open IMEM input file.";
		imem.close();
	}

	bitset<32> readInstr(bitset<32> ReadAddress)
	{
		// read instruction memory
		// return bitset<32> val
		// Fill in
		unsigned int addr = ReadAddress.to_ulong();
		unsigned long instruction = 0;

		// Combine 4 bytes: IMem[addr] is the most significant byte
		instruction = (IMem[addr].to_ulong() << 24) |
					  (IMem[addr + 1].to_ulong() << 16) |
					  (IMem[addr + 2].to_ulong() << 8) |
					  (IMem[addr + 3].to_ulong());

		return bitset<32>(instruction);
	}

private:
	vector<bitset<8>> IMem;
};

class DataMem
{
public:
	string id, opFilePath, ioDir;
	DataMem(string name, string ioDir) : id{name}, ioDir{ioDir}
	{
		DMem.resize(MemSize);
		opFilePath = ioDir + "\\" + name + "_DMEMResult.txt";
		ifstream dmem;
		string line;
		int i = 0;
		dmem.open(ioDir + "\\dmem.txt");
		if (dmem.is_open())
		{
			while (getline(dmem, line))
			{
				DMem[i] = bitset<8>(line);
				i++;
			}
		}
		else
			cout << "Unable to open DMEM input file.";
		dmem.close();
	}

	bitset<32> readDataMem(bitset<32> Address)
	{
		// read data memory
		// return bitset<32> val
		// Fill in
		unsigned int addr = Address.to_ulong();
		unsigned long data = 0;
		// Combine 4 bytes: DMem[addr] is the most significant byte
		data = (DMem[addr].to_ulong() << 24) |
			   (DMem[addr + 1].to_ulong() << 16) |
			   (DMem[addr + 2].to_ulong() << 8) |
			   (DMem[addr + 3].to_ulong());
		return bitset<32>(data);
	}

	void writeDataMem(bitset<32> Address, bitset<32> WriteData)
	{
		// write into memory
		// Fill in
		unsigned int addr = Address.to_ulong();
		// Split 4 bytes: DMem[addr] is the most significant byte
		DMem[addr] = bitset<8>(WriteData.to_ulong() >> 24);
		DMem[addr + 1] = bitset<8>(WriteData.to_ulong() >> 16);
		DMem[addr + 2] = bitset<8>(WriteData.to_ulong() >> 8);
		DMem[addr + 3] = bitset<8>(WriteData.to_ulong());
	}

	void outputDataMem()
	{
		ofstream dmemout;
		dmemout.open(opFilePath, std::ios_base::trunc);
		if (dmemout.is_open())
		{
			for (int j = 0; j < 1000; j++)
			{
				dmemout << DMem[j] << endl;
			}
		}
		else
			cout << "Unable to open " << id << " DMEM result file." << endl;
		dmemout.close();
	}

private:
	vector<bitset<8>> DMem;
};

class RegisterFile
{
public:
	string outputFile;
	RegisterFile(string ioDir) : outputFile{ioDir + "RFResult.txt"}
	{
		Registers.resize(32);
		Registers[0] = bitset<32>(0);
	}

	bitset<32> readRF(bitset<5> Reg_addr)
	{
		// Fill in
		unsigned int addr = Reg_addr.to_ulong();
		if (addr < 32)
		{
			return Registers[addr]; // Return the value of the register at addrss
		}
		else
		{
			cout << "Invalid register address: " << addr << endl;
			return bitset<32>(0); // Return 0 if the address is invalid
		}
	}

	void writeRF(bitset<5> Reg_addr, bitset<32> Wrt_reg_data)
	{
		// Fill in
		unsigned int addr = Reg_addr.to_ulong();
		if (addr > 0 && addr < 32) // Do not write to register 0
		{
			Registers[addr] = Wrt_reg_data; // Write the value to the register at address
		}
		else
		{
			cout << "Invalid register address: " << addr << endl;
			// Do nothing if the address is invalid
		}
	}

	void outputRF(int cycle)
	{
		ofstream rfout;
		if (cycle == 0)
			rfout.open(outputFile, std::ios_base::trunc); // create a new file
		else
			rfout.open(outputFile, std::ios_base::app); // append to the existing file
		if (rfout.is_open())
		{
			rfout << "----------------------------------------------------------------------" << endl;
			rfout << "State of RF after executing cycle:\t" << cycle << endl;
			for (int j = 0; j < 32; j++)
			{
				rfout << Registers[j] << endl;
			}
		}
		else
			cout << "Unable to open RF output file." << endl;
		rfout.close();
	}

private:
	vector<bitset<32>> Registers;
};


bitset<32> ALUOperation (short ALUOP, bitset<32> oprand1, bitset<32> oprand2)
{   
int result; 
	switch(ALUOP) 
		{							
		case add : result = oprand1.to_ulong() + oprand2.to_ulong(); break; // add
		case sub : result = oprand1.to_ulong() - oprand2.to_ulong(); break; // sub
		case and : result = oprand1.to_ulong() & oprand2.to_ulong(); break; // and
		case or  : result = oprand1.to_ulong() | oprand2.to_ulong(); break; // or
		case xor : result = oprand1.to_ulong() ^ oprand2.to_ulong(); break; // xor
		default : cout << "ALU operation not supported" << endl; break; // default case
		} 
		
	return bitset<32>(result);
}            

// Sign-extend a 12-bit input to 32 bits.
bitset<32> signExtend_32(const bitset<12>& input) {
    unsigned long value = input.to_ulong();
    // Check if the sign bit (bit 11) is 1.
    if (input.test(11)) {
        // Set the upper 20 bits to 1.
        value |= (~0UL << 12);
    }
    return bitset<32>(value);
}

// Sign-extend a 16-bit input to 32 bits.
bitset<32> signExtend_32(const bitset<16>& input) {
	unsigned long value = input.to_ulong();
	// Check if the sign bit (bit 15) is 1.
	if (input.test(15)) {
		// Set the upper 16 bits to 1.
		value |= (~0UL << 16);
	}
	return bitset<32>(value);
}


// Sign-extend a 20-bit input to 32 bits.
bitset<32> signExtend_32(const bitset<20>& input) {
    unsigned long value = input.to_ulong();
    // Check if the sign bit (bit 19) is 1.
    if (input.test(19)) {
        // Set the upper 12 bits to 1.
        value |= (~0UL << 20);
    }
    return bitset<32>(value);
}

bitset<16> signExtend_16(const bitset<12>& input) {
	unsigned long value = input.to_ulong();
	// Check if the sign bit (bit 11) is 1.
	if (input.test(11)) {
		// Set the upper 20 bits to 1.
		value |= (~0UL << 12);
	}
	return bitset<16>(value);
}

class Core
{
public:
	RegisterFile myRF;
	uint32_t cycle = 0;
	bool halted = false;
	string ioDir;
	struct stateStruct state, nextState;
	InsMem ext_imem;
	DataMem ext_dmem;

	Core(string ioDir, InsMem &imem, DataMem &dmem) : myRF(ioDir), ioDir{ioDir}, ext_imem{imem}, ext_dmem{dmem} {}

	virtual void step() {}

	virtual void printState() {}
};

class SingleStageCore : public Core
{
public:
	short  ALUOP = 0; // ALU operation code
	SingleStageCore(string ioDir, InsMem &imem, DataMem &dmem) : Core(ioDir + "\\SS_", imem, dmem), opFilePath(ioDir + "\\StateResult_SS.txt") {
		// Initialize nextState structures
		state.IF.PC = bitset<32>(0);
		state.IF.nop = false;
		state.IF.nop = false;
	}

	void step()
	{

		/* Your implementation*/
		// 1. IF: Fetch the instruction from the instruction memory using the current PC value.
		bitset<32> instruction = ext_imem.readInstr(state.IF.PC);
		if (instruction == 0xffffffff) // halt instruction
            {
				nextState.IF.PC = state.IF.PC;
				nextState.IF.nop = 1;

				if (state.IF.nop)
					halted = true;
				myRF.outputRF(cycle);		  // dump RF
				printState(nextState, cycle); // print states after executing cycle 0, cycle 1, cycle 2 ...

				state = nextState; // The end of the cycle and updates the current state with the values calculated in this cycle
				cycle++;
				return;                            
			}
		else
			{
				nextState.IF.PC = state.IF.PC.to_ulong() + 4; // Increment PC for the next instruction
				nextState.IF.nop = 0; // Valid instruction
			}

		// 2. ID: Decode the instruction and update the ID stage state.
		// Extract the opcode and other fields from the instruction
		string opcode = instruction.to_string().substr(25, 7); // last 7 bits
		bool is_R_type = false; // R-type instruction flag
		bool is_I_type = false; // I-type instruction flag
		// Read the register file to get the values of the source registers.
		// Identify the instruction type (R, I, S, B, J) based on the opcode, specifically, the instructions are:
		// 0110011: R-type (add, sub, and, or, xor)
		// 0010011: I-type (addi, andi, ori, xori)
		// 0000011: I-type (lw)
		// 0100011: S-type (sw)
		// 1100011: B-type (beq, bne)
		// 1101111: J-type (jal)

		string instrStr = instruction.to_string(); // Convert bitset to string for easier manipulation
		if (opcode == "0110011") { 
			// -------- R-type --------
			// R-type format (RISC-V): 
			// [31:25] funct7, [24:20] rs2, [19:15] rs1, [14:12] funct3, [11:7] rd, [6:0] opcode
			// Extract register fields:
			string rd_str  = instrStr.substr(20, 5); // bits 11:7
			string rs1_str = instrStr.substr(12, 5); // bits 19:15
			string rs2_str = instrStr.substr(7, 5);  // bits 24:20
			bitset<5> rd(rd_str), rs1(rs1_str), rs2(rs2_str);
			// Read source registers.
			bitset<32> val1 = myRF.readRF(rs1);
			bitset<32> val2 = myRF.readRF(rs2);
			
			// Update EX stage state.
			state.EX.Read_data1 = val1;
			state.EX.Read_data2 = val2;
			state.EX.nop = false;
			// Determine ALU operation based on funct3 and funct7
			// funct3 is bits 17:19, fuct7 is bits 0:6, big-endian
			bitset<3> funct3 = bitset<3>(instrStr.substr(17, 3)); // bits 14:12
			bitset<7> funct7 = bitset<7>(instrStr.substr(0, 7)); // bits 31:25
			if (funct3 == bitset<3>(0b000) && funct7 == bitset<7>(0b0000000)) {
				ALUOP = add; // add
			}
			else if (funct3 == bitset<3>(0b000) && funct7 == bitset<7>(0b0100000)) {
				ALUOP = sub; // sub
			}
			else if (funct3 == bitset<3>(0b111)) {
				ALUOP = and; // and
			}
			else if (funct3 == bitset<3>(0b110)) {
				ALUOP = or;  // or
			}
			else if (funct3 == bitset<3>(0b100)) {
				ALUOP = xor; // xor
			}
			else {
				cout << "Unsupported R-type instruction" << endl;
				nextState.IF.nop = true; // Set nop to true for unsupported instructions
			}
			state.EX.nop = false;
			state.MEM.nop = true;
			state.WB.nop = false;
			state.EX.Wrt_reg_addr = rd; // Set the destination register for write-back

			is_R_type = true; // Set R-type flag
		}
		else if (opcode == "0010011") { 
			// -------- I-type--------
			// I-type format: [31:20] immediate, [19:15] rs1, [14:12] funct3, [11:7] rd, [6:0] opcode
			string rd_str  = instrStr.substr(20, 5);  // bits 11:7
			string rs1_str = instrStr.substr(12, 5);  // bits 19:15
			string imm_str = instrStr.substr(0, 12);    // bits 31:20 (12-bit immediate)
			bitset<5> rd(rd_str), rs1(rs1_str);
			// For simplicity, treat immediate as unsigned.
			bitset<32> val1 = myRF.readRF(rs1);
			
			// Update EX stage state.
			state.EX.Read_data1 = val1;
			state.EX.Read_data2 = signExtend_32(bitset<12>(imm_str)); // Sign-extend the immediate value
			state.EX.nop = false;
			// ALU operation for I-type instructions (e.g., addi, andi, ori, xori)
			bitset<3> funct3 = bitset<3>(instrStr.substr(17, 3)); // bits 14:12
			if (funct3 == bitset<3>(0b000)) {
				ALUOP = add; // addi
			}
			else if (funct3 == bitset<3>(0b111)) {
				ALUOP = and; // andi
			}
			else if (funct3 == bitset<3>(0b110)) {
				ALUOP = or;  // ori
			}
			else if (funct3 == bitset<3>(0b100)) {
				ALUOP = xor; // xori
			}
			else {
				cout << "Unsupported I-type instruction" << endl;
				nextState.IF.nop = true; // Set nop to true for unsupported instructions
			}
			state.EX.nop = false;
			state.MEM.nop = true;
			state.WB.nop = false;
			state.EX.Wrt_reg_addr = rd; // Set the destination register for write-back

			is_I_type = true; // Set I-type flag
			
		}
		else if (opcode == "0000011") { 
			// -------- I-type: lw (load word) --------
			// Format same as I-type arithmetic.
			string rd_str  = instrStr.substr(20, 5);
			string rs1_str = instrStr.substr(12, 5);
			string imm_str = instrStr.substr(0, 12);
			bitset<5> rd(rd_str), rs1(rs1_str);
			
			// Update EX stage state.
			bitset<32> val1 = myRF.readRF(rs1);
			state.EX.Read_data1 = val1;
			state.EX.Read_data2 = signExtend_32(bitset<12>(imm_str)); // Sign-extend the immediate value
			state.EX.nop = false;

			// ALU operation for lw: add base address and immediate to get the effective address.
			ALUOP = add; // lw uses addition

			state.EX.nop = false;
			state.EX.Wrt_reg_addr = rd; // Write back to rd register
			// Set MEM stage state for lw.
			state.MEM.nop = false; // Not a nop instruction
			state.MEM.rd_mem = true; // Set read memory flag
			state.MEM.wrt_mem = false; // Not a write memory operation
			state.WB.nop = false; // Not a nop instruction
		}
		else if (opcode == "0100011") { 
			// -------- S-type: sw (store word) --------
			// S-type format: [31:25] immediate[11:5] | [11:7] immediate[4:0], [19:15] rs1, [24:20] rs2, [6:0] opcode
			string imm_high = instrStr.substr(0, 7);   // bits 31:25
			string imm_low  = instrStr.substr(20, 5);   // bits 11:7
			string rs1_str  = instrStr.substr(12, 5);    // bits 19:15
			string rs2_str  = instrStr.substr(7, 5);     // bits 24:20
			string imm_str  = imm_high + imm_low;         // 12-bit immediate
			bitset<5> rs1(rs1_str), rs2(rs2_str);

			// Update EX stage state.
			bitset<32> val1 = myRF.readRF(rs1);
			bitset<32> val2 = signExtend_32(bitset<12>(imm_str));
			state.EX.Read_data1 = val1;
			state.EX.Read_data2 = val2; // val2 is the immediate value
			state.EX.nop = false;
			// ALU operation for sw: add base address and immediate to get the effective address.
			ALUOP = add; // sw uses addition

			state.EX.nop = false;
			// Set MEM stage state for sw.
			state.MEM.nop = false; // Not a nop instruction
			state.MEM.rd_mem = false; // Not a read memory operation
			state.MEM.wrt_mem = true; // Set write memory flagstate.MEM.Wrt_reg_addr = rs2; // Write data from rs2 register
			state.MEM.Store_data = myRF.readRF(rs2); // Store data to be written to memory
			state.WB.nop = true; // No write-back for sw
		}
		else if (opcode == "1100011") { 
			// -------- B-type: beq (branch if equal) --------
			string rs1_str = instrStr.substr(12, 5); // bits 19:15
			string rs2_str = instrStr.substr(7, 5);  // bits 24:20
			// In a real implementation, the branch offset would be decoded from various bits.
			bitset<5> rs1(rs1_str), rs2(rs2_str);
			bitset<32> val1 = myRF.readRF(rs1);
			bitset<32> val2 = myRF.readRF(rs2);
			// immediate for B type:
			// imm[12]: Bit 31 of the instruction
			// imm[10:5]: Bits 30–25 of the instruction
			// imm[4:1]: Bits 11–8 of the instruction
			// imm[11]: Bit 7 of the instruction
			bitset<12> imm = bitset<12>(instrStr.substr(0, 1) + instrStr.substr(24, 1) + instrStr.substr(1, 6) + instrStr.substr(20, 4));
			int branchOffset = signExtend_32(imm).to_ulong() << 1; // Sign-extend the immediate value
			// beq and bne:
			// beq: branch if equal, bne: branch if not equal
			// funct3 is bits 14:12
			bitset<3> funct3 = bitset<3>(instrStr.substr(17, 3)); // bits 14:12
			if (funct3 == bitset<3>(0b000)) { // beq
				if (val1 == val2) {
					// Update PC with branch target address.
					nextState.IF.PC = bitset<32>(state.IF.PC.to_ulong() + branchOffset);
				}
			}
			else if (funct3 == bitset<3>(0b001)) { // bne
				if (val1 != val2) {
					// Update PC with branch target address.
					nextState.IF.PC = bitset<32>(state.IF.PC.to_ulong() + branchOffset);
				}
			}
			else {
				cout << "Unsupported B-type instruction" << endl;
				nextState.IF.nop = true; // Set nop to true for unsupported instructions
			}
			state.EX.nop = true; // No operation in EX stage for branch instructions
			state.MEM.nop = true; // No operation in MEM stage for branch instructions
			state.WB.nop = true; // No operation in WB stage for branch instructions
		}
		else if (opcode == "1101111") { 
			// -------- J-type: jal (jump and link) --------
			// J-type format: [31:12] immediate, [11:7] rd, [6:0] opcode
			string rd_str  = instrStr.substr(20, 5);   // bits 11:7
			// immediate for J type: 20, 10:1, 11, 19:12
			string imm_str = instrStr.substr(0, 1)       // imm[20]
				+ instrStr.substr(12, 8)     // imm[19:12]	
				+ instrStr.substr(11, 1)      // imm[11]
                + instrStr.substr(1, 10)      // imm[10:1]
            ;
			int branchOffset = signExtend_32(bitset<20>(imm_str)).to_ulong() << 1; // Sign-extend the immediate value
			bitset<5> rd(rd_str);
			// Save return address (current PC + 4) in rd.
			myRF.writeRF(rd, bitset<32>(state.IF.PC.to_ulong() + 4));
			// Update PC with jump target.
			nextState.IF.PC = bitset<32>(state.IF.PC.to_ulong() + branchOffset);
			state.EX.nop = true;
			state.MEM.nop = true;
			state.WB.nop = true;
		}
		else {
			cout << "Unsupported instruction" << endl;
			state.EX.nop = true;
			state.MEM.nop = true;
			state.WB.nop = true;
		}


		// 3. EX: Perform the ALU operation based on the instruction type.
		if (state.EX.nop == false) {
			// ALU operation
			state.MEM.ALUresult = ALUOperation(ALUOP, state.EX.Read_data1, state.EX.Read_data2);
			if(is_R_type || is_I_type) 
				state.WB.Wrt_data = state.MEM.ALUresult; // Set the write-back data to the ALU result
		}

		// 4. MEM: Access the data memory if needed.
		if (state.MEM.nop == false) {
			if (state.MEM.rd_mem) {
				// R and I type instructions
				state.WB.Wrt_data = ext_dmem.readDataMem(state.MEM.ALUresult);
			}
			else if (state.MEM.wrt_mem) {
				// Store word to memory
				ext_dmem.writeDataMem(state.MEM.ALUresult, state.MEM.Store_data);
			}
		}

		// 5. WB: Write back the result to the register file if needed.
		if (state.WB.nop == false) {
			state.WB.Wrt_reg_addr = state.EX.Wrt_reg_addr; // Set the destination register for write-back
			myRF.writeRF(state.WB.Wrt_reg_addr, state.WB.Wrt_data);
		}


		myRF.outputRF(cycle);		  // dump RF
		printState(nextState, cycle); // print states after executing cycle 0, cycle 1, cycle 2 ...

		state = nextState; // The end of the cycle and updates the current state with the values calculated in this cycle
		cycle++;
	}

	void printState(stateStruct state, int cycle)
	{
		ofstream printstate;
		if (cycle == 0)
			printstate.open(opFilePath, std::ios_base::trunc);
		else
			printstate.open(opFilePath, std::ios_base::app);
		if (printstate.is_open())
		{	
			printstate << "----------------------------------------------------------------------" << endl;		
			printstate << "State after executing cycle:\t" << cycle << endl;
			printstate << "IF.PC:\t" << state.IF.PC.to_ulong() << endl;
			printstate << "IF.nop:\t" << state.IF.nop << endl;
		}
		else
			cout << "Unable to open SS StateResult output file." << endl;
		printstate.close();
	}

private:
	string opFilePath;
};

class FiveStageCore : public Core
{
public:
	short  ALUOP = 0; // ALU operation code
	bitset<32> instruction = bitset<32>(0); // Initialize instruction to 0
	FiveStageCore(string ioDir, InsMem &imem, DataMem &dmem) : Core(ioDir + "\\FS_", imem, dmem), opFilePath(ioDir + "\\StateResult_FS.txt") {
		// Initialize state structures
		// IF
		state.IF.PC = bitset<32>(0);
		state.IF.nop = false;
		// ID
		state.ID.Instr = bitset<32>(0);
		state.ID.nop = true;
		// EX
		state.EX.Read_data1 = bitset<32>(0);
		state.EX.Read_data2 = bitset<32>(0);
		state.EX.Imm = bitset<16>(0);
		state.EX.Rs = bitset<5>(0);
		state.EX.Rt = bitset<5>(0);
		state.EX.Wrt_reg_addr = bitset<5>(0);
		state.EX.is_I_type = false;
		state.EX.alu_op = false;
		state.EX.rd_mem = false;
		state.EX.wrt_mem = false;
		state.EX.wrt_enable = false;
		state.EX.nop = true;
		// MEM
		state.MEM.ALUresult = bitset<32>(0);
		state.MEM.Store_data = bitset<32>(0);
		state.MEM.Rs = bitset<5>(0);
		state.MEM.Rt = bitset<5>(0);
		state.MEM.Wrt_reg_addr = bitset<5>(0);
		state.MEM.wrt_enable = false;
		state.MEM.rd_mem = false;
		state.MEM.wrt_mem = false;
		state.MEM.nop = true;
		// WB
		state.WB.Wrt_data = bitset<32>(0);
		state.WB.Wrt_reg_addr = bitset<5>(0);
		state.WB.wrt_enable = false;
		state.WB.nop = true;
	}

	void step()
	{
		nextState = state; // Initialize nextState to the current state+
		
		/* Your implementation */
		/* --------------------- WB stage --------------------- */
		// 1. WB: Write back the result to the register file if needed.
		if (state.WB.nop == false)
		{
			if (state.WB.wrt_enable){
				myRF.writeRF(state.WB.Wrt_reg_addr, state.WB.Wrt_data);
				nextState.WB.wrt_enable = state.MEM.wrt_enable; // Set write enable for WB stage based on MEM stage			
			}
			nextState.WB.nop = state.MEM.nop; // Set nop for WB stage based on MEM stage				
		} 

		/* --------------------- MEM stage -------------------- */
		// 1. MEM: Access the data memory if needed.
		// Memory stage: for a load, read data; for a store, write data; for add, just pass ALU result.
		if (state.MEM.nop == false)
		{
			if (state.MEM.rd_mem)  // lw instruction
			{
				nextState.WB.Wrt_data = ext_dmem.readDataMem(state.MEM.ALUresult);
				nextState.WB.wrt_enable = true;
				nextState.WB.Wrt_reg_addr = state.MEM.Wrt_reg_addr; // Set the destination register for write-back
			}
			else if (state.MEM.wrt_mem)  // sw instruction
			{
				ext_dmem.writeDataMem(state.MEM.ALUresult, state.MEM.Store_data);
				nextState.WB.wrt_enable = false;
			}
			else
			{
				nextState.WB.Wrt_data = state.MEM.ALUresult;	
				nextState.WB.wrt_enable = state.MEM.wrt_enable;
				nextState.WB.Wrt_reg_addr = state.MEM.Wrt_reg_addr; // Restore the line to set the destination register for write-back
			}
			// nextState.WB.Wrt_reg_addr = state.MEM.Wrt_reg_addr;
			// nextState.MEM.Wrt_reg_addr = state.EX.Wrt_reg_addr; // Set the destination register for write-back
			
			nextState.WB.nop = false; // Not a nop instruction
			nextState.MEM.nop = state.EX.nop; // Set nop for MEM stage based on ID stage
		}

		/* --------------------- EX stage --------------------- */
		// Execute stage: perform ALU operations (for add, compute sum; for ld/sw, compute effective address).
		if (state.EX.nop == false)
		{
			if (state.EX.is_I_type && state.EX.wrt_mem) { 
				// For sw instruction, compute effective address.
				int imm = signExtend_32(state.EX.Imm).to_ulong(); // Sign-extend the immediate value
				nextState.MEM.ALUresult = state.EX.Read_data1.to_ulong() + imm;
				nextState.MEM.Store_data = state.EX.Read_data2;

				// For sw, the data to store is in Read_data2.
				nextState.MEM.Wrt_reg_addr = state.EX.Wrt_reg_addr; // Set the destination register for write-back
				nextState.MEM.wrt_enable = state.EX.wrt_enable;
			}
			else if (state.EX.is_I_type && state.EX.rd_mem) { 
				// For lw instruction, compute effective address.
				int imm = signExtend_32(state.EX.Imm).to_ulong(); // Sign-extend the immediate value
				nextState.MEM.ALUresult = state.EX.Read_data1.to_ulong() + imm;
			}
			else if (state.EX.is_I_type) { // For I-type instructions (e.g., addi, andi, ori, xori)
				nextState.MEM.ALUresult = ALUOperation(ALUOP, state.EX.Read_data1, signExtend_32(state.EX.Imm)); // Sign-extend the immediate value
			}
			else { // For R-type instructions and J-type instructions
				nextState.MEM.ALUresult = ALUOperation(ALUOP, state.EX.Read_data1, state.EX.Read_data2);
			}
			nextState.MEM.Wrt_reg_addr = state.EX.Wrt_reg_addr; // Set the destination register for write-back
			nextState.MEM.wrt_enable = state.EX.wrt_enable;
			nextState.MEM.rd_mem = state.EX.rd_mem;
			nextState.MEM.wrt_mem = state.EX.wrt_mem;
			nextState.MEM.nop = false;
			nextState.MEM.nop = state.EX.nop; // Set nop for MEM stage based on EX stage
			nextState.MEM.Rs = state.EX.Rs; // Set Rs for MEM stage based on EX stage
			nextState.MEM.Rt = state.EX.Rt; // Set Rt for MEM stage based on EX stage
			nextState.EX.nop = (state.ID.nop && state.IF.nop); // Set nop for EX stage based on ID stage

		}

		

		/* --------------------- ID stage --------------------- */
		// 1. ID: Decode the instruction and update the ID stage state.
		// Extract the opcode and other fields from the instruction
		bool stall = false; // Initialize stall flag to false
		// Instruction Decode stage: decode the instruction from the IF/ID register.
		if (state.ID.nop == false)
		{
			string instrStr = state.ID.Instr.to_string();
			// Assume opcode is the last 7 bits (indices 25-31).
			string opcode = instrStr.substr(25, 7);

			// Detect Load-use hazard: if the instruction in EX stage is a load and the instruction in ID stage reads from the same register, stall the pipeline.
			if (state.EX.nop == false && state.EX.rd_mem)  // EX stage is a load
			{
				bitset<5> loadDest = state.EX.Wrt_reg_addr;
				if (opcode == "0110011" || opcode == "0100011") { // R and S type: source in bits 12-16 and 7-11.
					bitset<5> src1(instrStr.substr(12, 5));
					bitset<5> src2(instrStr.substr(7, 5));
					if (src1 == loadDest || src2 == loadDest)
						stall = true;
				}
				else if (opcode == "0000011") { // lw: source in bits 12-16.
					bitset<5> src(instrStr.substr(12, 5));
					if (src == loadDest)
						stall = true;
				}
				
			}
			if (stall)
			{
				// Load-use hazard detected: insert a stall.
				// Insert a bubble in the EX stage by setting nextState.EX.nop.
				nextState.EX.nop = true;
				// Do not advance the current instruction in ID; keep it in the same stage.
				// Also, stall the IF stage by keeping the PC the same.
				nextState.IF.PC = state.IF.PC; // Stall the IF stage.
				nextState.EX.is_I_type = false;
				nextState.EX.rd_mem = false;
				nextState.EX.wrt_mem = false;
				nextState.EX.wrt_enable = false;
				cout << "Load-use hazard detected: stalling pipeline." << endl;
			}else{
				// No hazard detected: proceed with instruction decode.
				nextState.EX.Wrt_reg_addr =  0;
				forwardingStruct fwd;

				bitset<5> src1, src2;
				if (opcode == "0110011" || opcode == "0100011" || opcode == "1100011") { // R and S and B type: source in bits 12-16 and 7-11.
					src1 = bitset<5>(instrStr.substr(12, 5));
					src2 = bitset<5>(instrStr.substr(7, 5));
				}
				else if (opcode == "0000011" || opcode == "0010011") { // lw and I type: source in bits 12-16.	
					src1 = bitset<5>(instrStr.substr(12, 5));
					// Only src1 is used for lw.
				}

				// Check forwarding from the EX stage (the instruction one stage ahead).
				if (state.EX.nop == false && state.EX.wrt_enable && state.EX.Wrt_reg_addr.to_ulong() != 0)
				{
					if (src1 == state.EX.Wrt_reg_addr)
						fwd.forwardA = 1; // Forward from EX stage.
					if ((opcode == "0110011" || opcode == "0100011" || opcode == "1100011") && src2 == state.EX.Wrt_reg_addr)
						fwd.forwardB = 1;
				}
				// Check forwarding from the MEM stage.
				if (state.MEM.nop == false && state.MEM.wrt_enable && state.MEM.Wrt_reg_addr.to_ulong() != 0)
				{
					if (src1 == state.MEM.Wrt_reg_addr && fwd.forwardA == 0)
						fwd.forwardA = 2; // Forward from MEM stage.
					if ((opcode == "0110011" || opcode == "0100011" || opcode == "1100011") && src2 == state.MEM.Wrt_reg_addr && fwd.forwardB == 0)
						fwd.forwardB = 2;
				}

				if (opcode == "0110011") { // R-type instruction
					nextState.EX.nop = false;
					nextState.EX.is_I_type = false;
					// For add, assume:
					//   rs in bits 12-16, rt in bits 17-21, rd in bits 22-26.
					nextState.EX.Rs = src1;
					nextState.EX.Rt = src2;
					nextState.EX.Wrt_reg_addr = bitset<5>(instrStr.substr(20, 5)); // bits 11-7
					nextState.EX.wrt_enable = true;
					nextState.EX.rd_mem = false;
					nextState.EX.wrt_mem = false;
					// Then apply forwarding if needed.
					bitset<32> regVal1 = myRF.readRF(src1);
					bitset<32> regVal2 = myRF.readRF(src2);
					if (fwd.forwardA == 1)
						regVal1 = nextState.MEM.ALUresult; // Forward from EX stage (via MEM stage).
					else if (fwd.forwardA == 2)
						regVal1 = nextState.WB.Wrt_data; // Forward from MEM/WB stage.
					if (fwd.forwardB == 1)
						regVal2 = nextState.MEM.ALUresult;
					else if (fwd.forwardB == 2)
						regVal2 = nextState.WB.Wrt_data;
						
					nextState.EX.Read_data1 = regVal1;
					nextState.EX.Read_data2 = regVal2;

					// Determine ALU operation based on funct3 and funct7
					// funct3 is bits 17:19, fuct7 is bits 0:6, big-endian
					bitset<3> funct3 = bitset<3>(instrStr.substr(17, 3)); // bits 14:12
					bitset<7> funct7 = bitset<7>(instrStr.substr(0, 7)); // bits 31:25
					if (funct3 == bitset<3>(0b000) && funct7 == bitset<7>(0b0000000)) {
						ALUOP = add; // add
					}
					else if (funct3 == bitset<3>(0b000) && funct7 == bitset<7>(0b0100000)) {
						ALUOP = sub; // sub
					}
					else if (funct3 == bitset<3>(0b111)) {
						ALUOP = and; // and
					}
					else if (funct3 == bitset<3>(0b110)) {
						ALUOP = or;  // or
					}
					else if (funct3 == bitset<3>(0b100)) {
						ALUOP = xor; // xor
					}
					else {
						cout << "Unsupported R-type instruction" << endl;
					}
					}
					// I-type:
					else if (opcode == "0010011") { // I instruction

						nextState.EX.nop = false;
						nextState.EX.is_I_type = true;
						// For addi, assume:
						//   rs (base) in bits 12-16, rd (destination) in bits 22-26, immediate in bits 0-11.
						nextState.EX.Rs = src1;
						nextState.EX.Wrt_reg_addr = bitset<5>(instrStr.substr(20, 5));
						nextState.EX.Imm = signExtend_16(bitset<12>(instrStr.substr(0, 12))); // Sign-extend the immediate value
						nextState.EX.wrt_enable = true;
						nextState.EX.rd_mem = false; // No memory read for addi
						nextState.EX.wrt_mem = false; // No memory write for addi
						// For addi, use forwarding for the base register if applicable.
						bitset<32> regVal = myRF.readRF(src1);

						if (fwd.forwardA == 1)
							regVal = nextState.MEM.ALUresult;
						else if (fwd.forwardA == 2)
							regVal = nextState.WB.Wrt_data;
						nextState.EX.Read_data1 = regVal;
						// Determine ALU operation based on funct3
						bitset<3> funct3 = bitset<3>(instrStr.substr(17, 3)); // bits 14:12
						if (funct3 == bitset<3>(0b000)) {
							ALUOP = add; // add
						}
						else if (funct3 == bitset<3>(0b111)) {
							ALUOP = and; // and
						}
						else if (funct3 == bitset<3>(0b110)) {
							ALUOP = or;  // or
						}
						else if (funct3 == bitset<3>(0b100)) {
							ALUOP = xor; // xor
						}
						else {
							cout << "Unsupported R-type instruction" << endl;
						}
					}
					else if (opcode == "0000011") { // lw instruction
						nextState.EX.nop = false;
						nextState.EX.is_I_type = true;
						// For ld, assume:
						//   rs (base) in bits 12-16, rd (destination) in bits 22-26, immediate in bits 0-11.
						nextState.EX.Rs = src1;
						nextState.EX.Wrt_reg_addr = bitset<5>(instrStr.substr(20, 5));
						nextState.EX.Imm = signExtend_16(bitset<12>(instrStr.substr(0, 12))); // Sign-extend the immediate value
						nextState.EX.wrt_enable = true;
						nextState.EX.rd_mem = true;  // lw reads from memory.
						nextState.EX.wrt_mem = false;
						// For lw, use forwarding for the base register if applicable.
						bitset<32> regVal = myRF.readRF(src1);
						if (fwd.forwardA == 1)
							regVal = nextState.MEM.ALUresult;
						else if (fwd.forwardA == 2)
							regVal = nextState.WB.Wrt_data;
						nextState.EX.Read_data1 = regVal;
						ALUOP = add; // lw uses addition
					}
					else if (opcode == "0100011") { // sw instruction
						nextState.EX.nop = false;
						nextState.EX.is_I_type = true;
						// For sw, assume:
						//   rs (base) in bits 12-16, rt (data source) in bits 17-21, immediate in bits 0-11.
						nextState.EX.Rs = src1;
						nextState.EX.Rt = src2;
						string imm_high = instrStr.substr(0, 7);   // bits 31:25
						string imm_low  = instrStr.substr(20, 5);   // bits 11:7
						nextState.EX.Imm = signExtend_16(bitset<12>(imm_high + imm_low));
						nextState.EX.wrt_enable = false;
						nextState.EX.rd_mem = false;
						nextState.EX.wrt_mem = true; // sw writes to memory.
						// For sw, apply forwarding for both base and data registers.
						bitset<32> baseVal = myRF.readRF(src1);
						bitset<32> dataVal = myRF.readRF(src2);
						if (fwd.forwardA == 1)
							baseVal = nextState.MEM.ALUresult;
						else if (fwd.forwardA == 2)
							baseVal = nextState.WB.Wrt_data;
						if (fwd.forwardB == 1)
							dataVal = nextState.MEM.ALUresult;
						else if (fwd.forwardB == 2)
							dataVal = nextState.WB.Wrt_data;
						nextState.EX.Read_data1 = baseVal;
						nextState.EX.Read_data2 = dataVal;
						ALUOP = add; // sw uses addition
					}
					// B-type: beq and bne
					else if (opcode == "1100011") { // beq and bne instructions
						nextState.EX.nop = false;
						// For beq/bne, assume:
						//   rs1 in bits 12-16, rs2 in bits 7-11, immediate in bits 0-11.
						bitset<12> imm = bitset<12>(instrStr.substr(0, 1) + instrStr.substr(24, 1) + instrStr.substr(1, 6) + instrStr.substr(20, 4));
						int branchOffset = signExtend_32(imm).to_ulong() << 1; // Sign-extend the immediate value
						bitset<32> regVal1 = myRF.readRF(src1);
						bitset<32> regVal2 = myRF.readRF(src2);
						if (fwd.forwardA == 1)
							regVal1 = state.MEM.ALUresult;
						else if (fwd.forwardA == 2)
							regVal1 = nextState.WB.Wrt_data;
						if (fwd.forwardB == 1)
							regVal2 = nextState.MEM.ALUresult;
						else if (fwd.forwardB == 2)
							regVal2 = nextState.WB.Wrt_data;
						string funct3 = instrStr.substr(17, 3);

						bool branchTaken = false;
						if (funct3 == "000") { // beq
							branchTaken = (regVal1.to_ulong() == regVal2.to_ulong());
						}
						else if (funct3 == "001") { // bne
							branchTaken = (regVal1.to_ulong() != regVal2.to_ulong());
						}
						else {
							cout << "Unsupported B-type instruction" << endl;
						}
						if (branchTaken) {
							// Update PC with branch target address.
							nextState.IF.PC = bitset<32>(state.IF.PC.to_ulong() + branchOffset - 4);
							stall = true; // Set stall flag to true for branch instructions
							nextState.ID.nop = true; // Set ID stage to nop
							nextState.EX.nop = true; // Set EX stage to nop
						}

					}
					// J-type: jal instruction
					else if (opcode == "1101111") { // jal instruction

						// J-type format: [31:12] immediate, [11:7] rd, [6:0] opcode
						string rd_str  = instrStr.substr(20, 5);   // bits 11:7
						// immediate for J type: 20, 10:1, 11, 19:12
						string imm_str = instrStr.substr(0, 1)       // imm[20]
							+ instrStr.substr(12, 8)     // imm[19:12]	
							+ instrStr.substr(11, 1)      // imm[11]
							+ instrStr.substr(1, 10)      // imm[10:1]
						;
						int branchOffset = signExtend_32(bitset<20>(imm_str)).to_ulong() << 1; // Sign-extend the immediate value
						bitset<5> rd(rd_str);
						// Offset
						nextState.EX.Read_data1 = bitset<32>(state.IF.PC.to_ulong()); // Save return address (current PC + 4) in rd.
						nextState.EX.Read_data2 = bitset<32>(0); // No second operand for jal
						nextState.EX.Wrt_reg_addr = rd; // Set the destination register for write-back
						
						// Update PC with jump target.
						nextState.IF.PC = bitset<32>(state.IF.PC.to_ulong() + branchOffset - 4);
						nextState.ID.nop = true; // Set ID stage to nop
						nextState.EX.is_I_type = false; // Set is_I_type to false for jal instruction
						stall = true; // Set stall flag to true for jal instruction


					}
					else {
						// Unknown opcode: treat as a no-op.
						nextState.ID.nop = true;
						cout << "Unknown opcode detected: stalling pipeline." << endl;
					}
			}
		}

		/* --------------------- IF stage --------------------- */
		// 1. IF: Fetch the instruction from the instruction memory using the current PC value.
		if(state.IF.nop == false)
		{			
			if(!stall){
				nextState.ID.nop = false; // Set IF stage to not a nop instruction
				unsigned long instrStr = ext_imem.readInstr(state.IF.PC).to_ulong(); // Read instruction from instruction memory
				if (instrStr == 0xffffffff) // halt instruction
				{
					nextState.IF.PC = state.IF.PC;
					nextState.IF.nop = 1;
					nextState.ID.nop = 1; // Set ID stage to nop
					nextState.EX.nop = state.EX.nop; // Set EX stage to nop
				}
				else
				{
					nextState.ID.Instr = bitset<32>(instrStr); // Update ID stage with the fetched instruction
					nextState.IF.PC = state.IF.PC.to_ulong() + 4; // Increment PC for the next instruction
					nextState.IF.nop = 0; // Valid instruction
				}
			} 	
		}
		

		if (state.IF.nop && state.ID.nop && state.EX.nop && state.MEM.nop && state.WB.nop)
			halted = true;

		myRF.outputRF(cycle);		  // dump RF
		printState(nextState, cycle); // print states after executing cycle 0, cycle 1, cycle 2 ...

		state = nextState; // The end of the cycle and updates the current state with the values calculated in this cycle
		cycle++;
	}

	void printState(stateStruct state, int cycle)
	{
		ofstream printstate;
		if (cycle == 0)
			printstate.open(opFilePath, std::ios_base::trunc);
		else
			printstate.open(opFilePath, std::ios_base::app);
		if (printstate.is_open())
		{
			printstate << "----------------------------------------------------------------------" << endl;
			printstate << "State after executing cycle:\t" << cycle << endl;

			printstate << "IF.nop:\t" << state.IF.nop << endl;
			printstate << "IF.PC:\t" << state.IF.PC.to_ulong() << endl;

			printstate << "ID.Instr:\t" << state.ID.Instr << endl;
			printstate << "ID.nop:\t" << state.ID.nop << endl;

			printstate << "EX.Read_data1:\t" << state.EX.Read_data1 << endl;
			printstate << "EX.Read_data2:\t" << state.EX.Read_data2 << endl;
			printstate << "EX.Imm:\t" << state.EX.Imm << endl;
			printstate << "EX.Rs:\t" << state.EX.Rs << endl;
			printstate << "EX.Rt:\t" << state.EX.Rt << endl;
			printstate << "EX.Wrt_reg_addr:\t" << state.EX.Wrt_reg_addr << endl;
			printstate << "EX.is_I_type:\t" << state.EX.is_I_type << endl;
			printstate << "EX.rd_mem:\t" << state.EX.rd_mem << endl;
			printstate << "EX.wrt_mem:\t" << state.EX.wrt_mem << endl;
			printstate << "EX.alu_op:\t" << state.EX.alu_op << endl;
			printstate << "EX.wrt_enable:\t" << state.EX.wrt_enable << endl;
			printstate << "EX.nop:\t" << state.EX.nop << endl;

			printstate << "MEM.ALUresult:\t" << state.MEM.ALUresult << endl;
			printstate << "MEM.Store_data:\t" << state.MEM.Store_data << endl;
			printstate << "MEM.Rs:\t" << state.MEM.Rs << endl;
			printstate << "MEM.Rt:\t" << state.MEM.Rt << endl;
			printstate << "MEM.Wrt_reg_addr:\t" << state.MEM.Wrt_reg_addr << endl;
			printstate << "MEM.rd_mem:\t" << state.MEM.rd_mem << endl;
			printstate << "MEM.wrt_mem:\t" << state.MEM.wrt_mem << endl;
			printstate << "MEM.wrt_enable:\t" << state.MEM.wrt_enable << endl;
			printstate << "MEM.nop:\t" << state.MEM.nop << endl;

			printstate << "WB.Wrt_data:\t" << state.WB.Wrt_data << endl;
			printstate << "WB.Rs:\t" << state.WB.Rs << endl;
			printstate << "WB.Rt:\t" << state.WB.Rt << endl;
			printstate << "WB.Wrt_reg_addr:\t" << state.WB.Wrt_reg_addr << endl;
			printstate << "WB.wrt_enable:\t" << state.WB.wrt_enable << endl;
			printstate << "WB.nop:\t" << state.WB.nop << endl;
		}
		else
			cout << "Unable to open FS StateResult output file." << endl;
		printstate.close();
	}

private:
	string opFilePath;
};

// Performance Matrix:
void printPerformanceMatrix(Core &ss, Core &fs, string ioDir)
{
	int ssCycle = ss.cycle;
	int fsCycle = fs.cycle;
	int ssNumInst = ss.ext_imem.size;
	int fsNumInst = fs.ext_imem.size;
	// Print performance matrix to PerformanceMatrix.txt
	ofstream perfMatrix;
	perfMatrix.open(ioDir + "\\PerformanceMatrix.txt", std::ios_base::trunc);
	if (perfMatrix.is_open())
	{
		perfMatrix << "Performance Matrix" << endl;
		perfMatrix << "----------------------------------------" << endl;
		perfMatrix << "Single Stage Core:" << endl;
		perfMatrix << "Cycles: " << ssCycle << endl;
		perfMatrix << "Instructions: " << ssNumInst << endl;
		perfMatrix << "CPI: " << (float)ssCycle / ssNumInst << endl;
		// Also IPC
		perfMatrix << "IPC: " << (float)ssNumInst / ssCycle << endl;
		perfMatrix << "----------------------------------------" << endl;
		perfMatrix << "Five Stage Core:" << endl;
		perfMatrix << "Cycles: " << fsCycle << endl;
		perfMatrix << "Instructions: " << fsNumInst << endl;
		perfMatrix << "CPI: " << (float)fsCycle / fsNumInst << endl;
		// Also IPC
		perfMatrix << "IPC: " << (float)fsNumInst / fsCycle << endl;
	}
	else
		cout << "Unable to open Performance Matrix output file." << endl;
	
}

int main(int argc, char *argv[])
{

	string ioDir = "";
	if (argc == 1)
	{
		cout << "Enter path containing the memory files: ";
		cin >> ioDir;
	}
	else if (argc > 2)
	{
		cout << "Invalid number of arguments. Machine stopped." << endl;
		return -1;
	}
	else
	{
		ioDir = argv[1];
		cout << "IO Directory: " << ioDir << endl;
	}

	InsMem imem = InsMem("Imem", ioDir);
	DataMem dmem_ss = DataMem("SS", ioDir);
	DataMem dmem_fs = DataMem("FS", ioDir);

	SingleStageCore SSCore(ioDir, imem, dmem_ss);
	FiveStageCore FSCore(ioDir, imem, dmem_fs);

	while (1)
	{
		if (!SSCore.halted)
			SSCore.step();
		
		if (!FSCore.halted)
			FSCore.step();

		if (SSCore.halted && FSCore.halted)
			break;
		// if (FSCore.halted)
		// 	break;
	}

	// dump SS and FS data mem.
	SSCore.ext_dmem.outputDataMem();
	FSCore.ext_dmem.outputDataMem();
	printPerformanceMatrix(SSCore, FSCore, ioDir);

	return 0;
}