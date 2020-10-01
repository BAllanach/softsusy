CALCULATER IN C++


#include<iostream>
using namespace std;
int main(int argc, char const *argv[])
{ int a,b;
char s;
cout<<"tipe your first number"<<endl;cin>>a;
cout<<"tipe your second number"<<endl;cin>>b;
cout<<"tipewhat you want to do/,*,-,+ or you want to get the remainder for this tipe %"<<endl;cin>>s;

if(s=='+'){
    cout<<"your sum is"<<a+b<<endl;}
else if(s=='-')
{cout<<"your subtracted number  is"<<a-b<<endl;}  

else if(s=='/'){
    cout<<"your subtracted number  is   "<<a/b<<endl;}
else if(s=='*'){
    cout<<"your subtracted number  is   "<<a*b<<endl;}   
else if(s=='%'){
    cout<<"your subtracted number  is   "<<a%b<<endl;} 
    return 0;
}

