#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class InsMem
{
public:
	string id, ioDir;
	int size;
	InsMem(string name, string ioDir)
	{
		id = name;
		ifstream imem;
		string line;
		size = 0;
		imem.open(ioDir + "\\imem.txt");
		if (imem.is_open())
		{
			while (getline(imem, line))
			{
				size++;
			}
		}
		else
			cout << "Unable to open IMEM input file.";
		imem.close();
	}
};

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
    // Initialize the instruction memory with the given directory
	InsMem imem = InsMem("Imem", ioDir);
}