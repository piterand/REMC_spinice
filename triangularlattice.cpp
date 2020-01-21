#include "triangularlattice.h"

Triangularlattice::Triangularlattice()
{

}

void Triangularlattice::dropLattice(double l, Vect size)
{
    this->clear();
    totalSize=size*l;
    double p[][2]={
        {l/2., 0},
        {l/4., l*sqrt(3)/4.},
        {l*3/4., l*sqrt(3)/4.}
    };
    double m[][2]={
        {1,0},
        {1/2., sqrt(3)/2.},
        {1/2., -sqrt(3)/2.}
    };

    //vector<Part*> top(size.x*4),bottom(size.x*4),left(size.y*4),right(size.y*4),front(size.z*4),back(size.z*4);
    double lx,ly;

    for (int i=0; i<size.x; i++){
        for (int j=0; j<size.y; j++){
                lx=l*(double)i;
                ly=l*(double)j;
                for (int c=0;c<3;c++){
                    Part *temp = new Part();
                    temp->pos=Vect(p[c][0]+lx+(ly*1/2.),p[c][1]+(ly*sqrt(3)/2.),0);
                    temp->m=Vect(m[c][0],m[c][1],0);
                    this->add(temp);
                }
                //верхняя грань
                if(i==size.x-1){

                        Part *temp = new Part();
                        temp->pos=Vect(p[1][0]+lx+l+(ly*1/2.),p[1][1]+(ly*sqrt(3)/2.),0);
                        temp->m=Vect(m[1][0],m[1][1],0);
                        this->add(temp);


                }

                if(j==size.y-1){

                        Part *temp = new Part();
                        temp->pos=Vect(p[0][0]+lx+(ly*1/2.)+l*1/2.,p[0][1]+(ly*sqrt(3)/2.)+l*sqrt(3)/2.,0);
                        temp->m=Vect(m[0][0],m[0][1],0);
                        this->add(temp);


                }


        }
    }

    this->setInteractionRange(0); //this->setInteractionRange(l/sqrt(8)*1.05);
}

void Triangularlattice::dropLattice2(double l, Vect size)
{
    this->clear();
    totalSize=size*l;
    double p[][2]={
        {l/2., 0},
        {l/4., l*sqrt(3)/4.},
        {l*3/4., l*sqrt(3)/4.}
    };
    double m[][2]={
        {1,0},
        {1/2., sqrt(3)/2.},
        {1/2., -sqrt(3)/2.}
    };

    //vector<Part*> top(size.x*4),bottom(size.x*4),left(size.y*4),right(size.y*4),front(size.z*4),back(size.z*4);
    double lx,ly;

    for (int i=0; i<size.x; i++){
        for (int j=0; j<size.y; j++){
                lx=l*(double)i;
                ly=l*(double)j;
                if(j%2==0){
                    for (int c=0;c<3;c++){
                        Part *temp = new Part();
                        temp->pos=Vect(p[c][0]+lx+l*1/2.,p[c][1]+(ly*sqrt(3)/2.),0);
                        temp->m=Vect(m[c][0],m[c][1],0);
                        this->add(temp);
                    }
                }
                else{
                    for (int c=0;c<3;c++){
                        Part *temp = new Part();
                        temp->pos=Vect(p[c][0]+lx,p[c][1]+(ly*sqrt(3)/2.),0);
                        temp->m=Vect(m[c][0],m[c][1],0);
                        this->add(temp);
                    }
                }

                //верхняя грань
//                if(i==size.x-1){

//                        Part *temp = new Part();
//                        temp->pos=Vect(p[1][0]+lx+l+(ly*1/2.),p[1][1]+(ly*sqrt(3)/2.),0);
//                        temp->m=Vect(m[1][0],m[1][1],0);
//                        this->add(temp);


//                }

//                if(j==size.y-1){

//                        Part *temp = new Part();
//                        temp->pos=Vect(p[0][0]+lx+(ly*1/2.)+l*1/2.,p[0][1]+(ly*sqrt(3)/2.)+l*sqrt(3)/2.,0);
//                        temp->m=Vect(m[0][0],m[0][1],0);
//                        this->add(temp);


//                }


        }
    }

    this->setInteractionRange(0); //this->setInteractionRange(l/sqrt(8)*1.05);
}

