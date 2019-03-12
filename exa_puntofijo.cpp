#include <iostream>
#include <iomanip>
#include <cmath>

double ref1(double x, double y){
    return std::sqrt( (y*y +1) /3);
}
double ref2(double x, double y){
    return std::sqrt( (x*x*x +5) /3*x);
}



int main(){
    double x{-2},y{2},temx{0},temy{0},
                dif1,dif2;
    int iter{13};
    /*
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
