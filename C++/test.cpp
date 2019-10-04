
#include <string>
#include <iostream>

using namespace std;

struct prova{
    int cas1;
    int cas2;
    int cas3;
};

int main(){
    prova prova1;
    prova1.cas1 = 2;
    string member = "cas1";
    cout << prova1.member << "\n";
}
