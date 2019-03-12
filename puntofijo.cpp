#include <iostream>
#include <iomanip>
#include <cmath>

double fun1(double x,double y){return (-x*x + 4*x + y - 0.5)/2.0;}
double fun2(double x,double y){
    return (-x*x - 4*y*y + 11*y + 4)/11.0;
}
double eje1(double x,double y){
    return sqrt(y*y-x+1);
}
double eje2(double x,double y){
    return (3*y+sin(x*x))/4.0;
}
double tar1(double x,double y){
    return (-x*x +y +2 +4*x)/4.0;
}
double tar2(double x,double y){
    return (-2*x*y +3 +5*y)/5.0;
}
double ref1(double x, double y){
    return (x*x-y*y-x-3+3*x)/3.0;
}
double ref2(double x, double y){
    return (y-x-1)/3.0;
}



int main(){
    double x{0},y{0},temx{0},temy{0},
                dif1,dif2;
    int iter{13};
    //*
    std::cout<<"Iter. simultaneas: "
                    <<std::endl;
    std::cout<<std::setprecision(6)
                << std::setw(3)
                << 0 <<":"
                << std::setw(9)
                << x <<" |"
                << std::setw(9)
                << y <<" | "<< "["
                << std::setw(9)
                << " "<<","
                << std::setw(9)
                << " " << " ]" << " |"
                << std::setw(9)
                << " "
                << std::endl;
    for(int i{1}; i<=iter; i++){
        temx=ref1(x,y);
        temy=ref2(x,y);
        dif1=temx -x;
        dif2=temy -y;
        x=temx;
        y=temy;
        std::cout<<std::setprecision(6)
                    << std::fixed
                    << std::setw(3)
                    << i <<":"
                    << std::setw(9)
                    << x<<" |"
                    << std::setw(9)
                    << y <<" | "<< "["
                    << std::setw(9)
                    << dif1<<","
                    << std::setw(9)
                    << dif2 << " ]" << " |"
                    << std::setw(9)
                    << std::abs( (std::abs(dif1)
                    < std::abs(dif2)) ? dif2 : dif1)
                    << std::endl;
    }
    /*/
    std::cout<<"Iter. sucesivas: "
                    <<std::endl;
    std::cout<<std::setprecision(6)
                << std::setw(3)
                << 0 <<":"
                << std::setw(9)
                << x <<" |"
                << std::setw(9)
                << y <<" | "<< "["
                << std::setw(9)
                << " "<<","
                << std::setw(9)
                << " " << " ]" << " |"
                << std::setw(9)
                << " "
                << std::endl;
    for(int i{1}; i<=iter; i++){
        temx=ref1(x,y);
        temy=ref2(temx,y);
        dif1=temx -x;
        dif2=temy -y;
        x=temx;
        y=temy;
        std::cout<<std::setprecision(6)
                    << std::fixed
                    << std::setw(3)
                    << i <<":"
                    << std::setw(9)
                    << x<<" |"
                    << std::setw(9)
                    << y <<" | "<< "["
                    << std::setw(9)
                    << dif1<<","
                    << std::setw(9)
                    << dif2 << " ]" << " |"
                    << std::setw(9)
                    << std::abs( (std::abs(dif1)
                    < std::abs(dif2)) ? dif2 : dif1)
                    << std::endl;
    }
    //*/

return 0;}
