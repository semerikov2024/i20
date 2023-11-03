#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
	if(argc==1)
	{
		cout<<"Usage:\n\t"<<argv[0]<<" filename\n";
		return 1;
	}

	ifstream f(argv[1], ios::in|ios::binary);
	if(!f)
	{
		cout<<"Can't open file "<<argv[1]<<"\n";
		return 2;
	}
	float number;
	while(f.read((char*)&number, 4))
		cout<<number<<endl;
}

