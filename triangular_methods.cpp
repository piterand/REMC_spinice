/*
 * triangular_methods.cpp
 *
 *          Replica MC modeling Triangular functions
 *
 * programmed by:
 *
 *          Andriushchenko Petr 2017.10.19
*/


void Triangular_create_center_of_circle(vector<double>& center_of_circle,PartArray &sys2,double &l){
    //1-3x+3x^2 количество гексагонов от длины стороны
    //cout<<sys2[0]->pos.x<<" "<<sys2[0]->pos.y<<endl;
///// ahtung, govnokod mode on //// тут я нахожу вернхюю левую частицу, и по ней определяю вертикальную левую частицу в первом ряду.
    int VL_particle=0;              // от нее дальше и идет построение.
    int blizh_sosed=0;
    double p1=1000000;
    double p2=1000000;
    double p3=1000000;
    int ll=0;               // размер цветочка=кол-во гексагонов на стороне

    std::vector<double> tempvec1;   // находим левый верхний спин
    for(unsigned i=0;i<sys2.size();i++){
        if(sys2[i]->pos.y<p1)
            p1=sys2[i]->pos.y;
    }
    for(unsigned i=0;i<sys2.size();i++){
        if(p1==sys2[i]->pos.y)
            tempvec1.push_back(i);
    }
    for(unsigned i=0;i<tempvec1.size();i++){
        if(sys2[tempvec1[i]]->pos.x<p2){
            p2=sys2[tempvec1[i]]->pos.x;
            VL_particle=tempvec1[i];
        }
    }
    // определение ближайшего соседа и определение l
    for(unsigned i=0;i<tempvec1.size();i++){
        if(abs(sys2[VL_particle]->pos.x-sys2[tempvec1[i]]->pos.x)<p3 && VL_particle!=tempvec1[i]){
            blizh_sosed=sys2[tempvec1[i]]->Id();
            p3=abs(sys2[VL_particle]->pos.x-sys2[tempvec1[i]]->pos.x);
        }
    }

    l=abs(sys2[VL_particle]->pos.x-sys2[blizh_sosed]->pos.x);
    l=round(l*10)/10.;  // определение l с точностью до 1 знака после запятой


//    Part *temp;
//    Vect tvect(sys2[VL_particle]->pos.x-sqrt(3.)*l*0.25,sys2[VL_particle]->pos.y+l*0.75,0);
//    temp = sys2.findByPosition(tvect,0.001);
//    VL_particle=temp->Id();


    ///// ahtung govnokod mode off bolee-menee ////
    ll=tempvec1.size();       // размер цветочка

    for(int vert=0;vert<2*ll;vert++){
        if(vert<ll){
            for(int horiz=0;horiz<ll+1+vert;horiz++){
                center_of_circle.push_back((sys2[VL_particle]->pos.x - 0.5 * l * (vert+1))+l*horiz);
                center_of_circle.push_back(sys2[VL_particle]->pos.y + sqrt(3.)/3.*l + l * sqrt(3.)/2. * (vert));
                if(horiz<ll+vert)
                {
                    center_of_circle.push_back((sys2[VL_particle]->pos.x - 0.5*l*(vert)) + l*horiz);
                    center_of_circle.push_back(sys2[VL_particle]->pos.y + sqrt(3.)/6.*l + l * sqrt(3.)/2. * (vert));
                }
            }
        }
        else
        {
            for(int horiz=0;horiz<ll+2*l+(2*l-vert);horiz++){
                if(horiz==0){}
                else{
                    center_of_circle.push_back((sys2[VL_particle]->pos.x - 0.5 * l * (2*l-vert-l)) - 2*l +l*horiz);
                    center_of_circle.push_back(sys2[VL_particle]->pos.y + sqrt(3.)/3.*l + l * sqrt(3.)/2. * (vert));
                }
                center_of_circle.push_back((sys2[VL_particle]->pos.x - 0.5 *l*(2*l-vert-l)) - 1.5*l + l*horiz);
                center_of_circle.push_back(sys2[VL_particle]->pos.y + sqrt(3.)/6.*l + l * sqrt(3.)/2. * (vert));

            }


        }
    }


}

