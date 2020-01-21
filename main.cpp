/*
 * Replica MC modeling HC functions
 *
 *
 * programmed by:
 *
 *          Andriushchenko Petr 2017.10.19
*/

#include <iostream>
#include <cmath>
#include <map>
#include "random.h"
#include "PartArray.h"
#include "config.h"
#include "squareisinglattice.h"

#include <random>
#include <iomanip>

#include <mpi.h>

//int seed = (time(0));
//int seed = 1;
//default_random_engine generator(seed);

double set_temp(int ,unsigned,double,double,int,double,double);
void HC_create_center_of_circle(vector<double>&,PartArray&,double&); // создание центров хонеакомбов
void Triangular_create_center_of_circle(vector<double>&,PartArray&,double&); // создание центров треугольной
bool form_cenetr_to_particles(double ,double , double ,double,double);       // проверка - являются ли соседями центры хонеакомбов
bool is_Next_for_cicle(double ,double , double ,double );            //
bool g_isNeighbours_for_vertex(double,double, double,double,double);
bool check_skalar(double,double, double ,double);
void ordering_vertex_spins(int, vector<double>& ,int ,PartArray &,double, unsigned);
int check_closed_vertex(int,int ,int*,PartArray &,unsigned);
int MaxClass(int , int* ,int* ,int);
int init_calc_max_cluster(int , int* ,int );
void MC_step(PartArray &, unsigned long long , double &,double );
void exact_dos(int*, int , int , int , PartArray );


double l; // параметр решетки вычисляется автоматически для гексагональной решетки

/////////////////////////
//char examplename_mfsys[50] = "lattice.mfsys"; // название загружаемой системы
unsigned type_of_system = 3; // !!!!!!!!!!!!!!!!!0=HC, 1=Tr, 2=Kag, 3=Square.
/////////////////////////

unsigned temperature_distribution = 0; // Распределение температур между репликами
// 0 - равномерно для линейного масштаба
// 1 - рвнометрно для логорифмического масштаба
// 2 - загрузка из файла

double mintemp=0.1;						// выбор min температуры
double maxtemp=5;						// выбор max температуры


bool replica_on=1;  //включить/выключить репличный обмен
bool out_current_values=0;  //выводить каждое значение

unsigned long long Prohod_MC_equilibration=1e4;            // МК прогрев
unsigned long long between_exchange_equilibration=1e3;     //между репличными обменами
unsigned long long Prohod_MC_sampling=1e7;                 // МК шагов сэмплирование
unsigned long long between_exchange=1e4;                   // Между обменами в сэмплировании
unsigned long long out_iterations=1e6;                     // количество вывода
unsigned long long iterations_wo_out=1;                    // МК внутри цикла mc_prohod
unsigned long long TOTAL_MC_Steps=iterations_wo_out*out_iterations;  // всего MK шагов

#include "honeycomb_methods.cpp"
#include "triangular_methods.cpp"

double set_temp(int i, unsigned temperature_distribution,double maxtemp, double mintemp, int size, double log_mintemp, double d_log_t){
    if(temperature_distribution==0)
        return ((i+1)*((maxtemp-mintemp)/size)+mintemp);
    else if (temperature_distribution==1)
        return (exp(log_mintemp+(d_log_t)*(i)));
    else if (temperature_distribution==2){
        char buff[50]; // буфер промежуточного хранения считываемого из файла текста

        ifstream fin("tm.dat"); // открыли файл для чтения
        for (int ii=0;ii<(4*i)+1;ii++){
            buff[0]=0;
            fin >> buff; // считали первое слово из файла

            ii++;
        }
        fin.close(); // закрываем файл
        return (double)atof(buff);
        //cout <<"!!!!!!!!!!!!!my rank" <<rank<<"     "<<temperature << endl;

    }
}

bool form_cenetr_to_particles(double x,double y, double x1,double y1,double center_to_particle_distance)
{
    return sqrt(pow((x1-x),2)+pow((y1-y),2))<(center_to_particle_distance*l*1.05);
}

bool is_Next_for_cicle(double x,double y, double x1,double y1)
{
    if(x==x1&&y==y1)
        return false;
    else
        return sqrt(pow((x1-x),2)+pow((y1-y),2))<((l*sqrt(3)/2.)*1.05);
}

bool g_isNeighbours_for_vertex(double x,double y, double x1,double y1,double center_to_center_distance)
{
    if(x==x1&&y==y1)
        return false;
    else
        return sqrt(pow((x1-x),2)+pow((y1-y),2))<(center_to_center_distance*l*1.05);
}

bool check_skalar(double x,double y, double x1,double y1)
{
    if((x*x1+y*y1)>0)
        return true;
    else
        return false;
}

void ordering_vertex_spins(int amount_of_circle, vector<double>& center_of_circle,int neighbours_for_cicle[][6],PartArray &sys2,double center_to_particle_distance, unsigned type_of_system){

    if(type_of_system==0){
        for(int i=0;i<amount_of_circle;i++){    // заполнение массива neighbours_for_cicle

            int jk=0;
            for(unsigned int j=0;j<sys2.size();j++){

                if (form_cenetr_to_particles(center_of_circle[2*i],center_of_circle[2*i+1],sys2[j]->pos.x,sys2[j]->pos.y,center_to_particle_distance)){
                    neighbours_for_cicle[i][jk]=j;
                    jk++;
                }
                if(jk==6)
                    break;

            }
        }
    }
    else if(type_of_system==1){
        for(int i=0;i<amount_of_circle;i++){    // заполнение массива neighbours_for_cicle

            int jk=0;
            for(unsigned int j=0;j<sys2.size();j++){

                if (form_cenetr_to_particles(center_of_circle[2*i],center_of_circle[2*i+1],sys2[j]->pos.x,sys2[j]->pos.y,center_to_particle_distance)){
                    neighbours_for_cicle[i][jk]=j;
                    jk++;
                }
                if(jk==3)
                    break;

            }
        }
    }

    //    if(type_of_system==0){
    //        cout<<"isNeighbours_for_cicle"<<endl;
    //        for(int ii=0;ii<amount_of_circle;ii++){
    //            for(int jj=0;jj<6;jj++){
    //                cout<<neighbours_for_cicle[ii][jj]<<",";
    //            }
    //            cout<<endl;
    //        }
    //    }
    //    if(type_of_system==1){
    //        cout<<"isNeighbours_for_cicle"<<endl;
    //        for(int ii=0;ii<amount_of_circle;ii++){
    //            for(int jj=0;jj<3;jj++){
    //                cout<<neighbours_for_cicle[ii][jj]<<",";
    //            }
    //            cout<<endl;
    //        }
    //    }

    // выравнивание по кругу
    if(type_of_system!=1)
    {
        int tmp;
        int t2;
        for(int i=0;i<amount_of_circle;i++){
            tmp=1;
            for(int j=0;j<6;j++){
                if(is_Next_for_cicle(sys2[neighbours_for_cicle[i][tmp-1]]->pos.x,sys2[neighbours_for_cicle[i][tmp-1]]->pos.y,sys2[neighbours_for_cicle[i][j]]->pos.x,sys2[neighbours_for_cicle[i][j]]->pos.y)){
                    t2=neighbours_for_cicle[i][tmp];
                    neighbours_for_cicle[i][tmp]=neighbours_for_cicle[i][j];
                    neighbours_for_cicle[i][j]=t2;
                    tmp++;
                    break;
                }
            }

            for(int j=tmp;j<6;j++){
                for(int k=tmp;k<6;k++){
                    if(is_Next_for_cicle(sys2[neighbours_for_cicle[i][tmp-1]]->pos.x,sys2[neighbours_for_cicle[i][tmp-1]]->pos.y,sys2[neighbours_for_cicle[i][k]]->pos.x,sys2[neighbours_for_cicle[i][k]]->pos.y)){
                        t2=neighbours_for_cicle[i][tmp];
                        neighbours_for_cicle[i][tmp]=neighbours_for_cicle[i][k];
                        neighbours_for_cicle[i][k]=t2;
                        tmp++;
                        break;


                    }

                }
            }


        }
    }

}

int check_closed_vertex(int amount_of_circle,int neighbours_for_cicle[][6],int* is_vertex_closed,PartArray &sys2,unsigned type_of_system){
    int amount=0;
    if(type_of_system==0){
        for(int i=0;i<amount_of_circle;i++){
            for(int j=0;j<6;j++){
                if(!check_skalar(sys2[neighbours_for_cicle[i][j]]->m.x,sys2[neighbours_for_cicle[i][j]]->m.y,sys2[neighbours_for_cicle[i][(int)fmod(j+1,6)]]->m.x,sys2[neighbours_for_cicle[i][(int)fmod(j+1,6)]]->m.y)){
                    is_vertex_closed[i]=0;
                    break;
                }

                if(j==5){
                    is_vertex_closed[i]=1;
                    amount++;
                }
            }
        }
    }
    else if(type_of_system==1){
        for(int i=0;i<amount_of_circle;i++){
            for(int j=0;j<3;j++){
                if(check_skalar(sys2[neighbours_for_cicle[i][j]]->m.x,sys2[neighbours_for_cicle[i][j]]->m.y,sys2[neighbours_for_cicle[i][(int)fmod(j+1,3)]]->m.x,sys2[neighbours_for_cicle[i][(int)fmod(j+1,3)]]->m.y)){
                    is_vertex_closed[i]=0;
                    break;
                }

                if(j==2){
                    is_vertex_closed[i]=1;
                    amount++;
                }
            }
        }
    }
    return amount;
}

int MaxClass(int per, int* counted_vertex,int* Ochered,int neighbours_for_vertex[][6]){
    int top;
    int w = 0;
    int r = 1;
    Ochered[0]=per;
    counted_vertex[per]=0;
    while(w<r){
        top=Ochered[w];

        for(int i=0;i<6;++i){
            if(counted_vertex[neighbours_for_vertex[top][i]]==1){
                Ochered[r]=neighbours_for_vertex[top][i];
                counted_vertex[neighbours_for_vertex[top][i]]=0;
                r++;
            }
            if(counted_vertex[neighbours_for_vertex[top][i]]==-1)
                break;

        }
        w++;
    }
    delete[] Ochered;
    return r;
}

int init_calc_max_cluster(int amount_of_circle, int* is_vertex_closed,int neighbours_for_vertex[][6] ){

    int Ochered[amount_of_circle];
    int Max_cluster=0;
    int t=0;
    for(int i=0;i<amount_of_circle;++i){
        if(is_vertex_closed[i]==1){
            int ttp = i;
            t = MaxClass(ttp, is_vertex_closed, Ochered, neighbours_for_vertex);
            if(t>Max_cluster)
                Max_cluster = t;

            t=0;


        }

    }
    return Max_cluster;
}

void MC_step(PartArray &sys2, unsigned long long iterations_wo_out, double &sysEnergy,double temperature){
    //uniform_int_distribution<int> distr_int(0,sys2.size()-1);
    //uniform_real_distribution<double> distr_double(0.,1.);

    double oldE,newE,delta_t;
    int num;

    for (unsigned long long step=0; step < iterations_wo_out; step++){
        //num = distr_int(generator);
        num=rand()%(sys2.size()-1);
        oldE = sysEnergy;

        delta_t=sys2.Check_dT(sys2[num]);
        newE = oldE+delta_t;
        double rand1=(double)(rand())/RAND_MAX;
        double veroyatnost=(exp(-(delta_t)/temperature));

        if (newE > oldE){
            if (veroyatnost<rand1){
                sysEnergy = oldE;
            }
            else {
                sysEnergy = newE;
                sys2[num]->rotate(true);
            }
        }
        else {
            sysEnergy = newE;
            sys2[num]->rotate(true);
        }
    }
}

void exact_dos(int* is_vertex_closed, int neighbours_for_cicle[][6], int neighbours_for_vertex[][6], int amount_of_circle, PartArray sys2)
{
    long e1=0;
    double cl_v=0;
    double m_cl=0;
    long long ti=0;

    sys2.state.reset();
    ofstream fq("dos_cl.csv");
    ofstream fqmcl("dos_mcl.csv");


    //struct dos3 {long e;int p};
    std::map<pair<long, int>, long > dos2_cl;
    std::map<pair<long, int>, long > dos2_mcl;

    //    dos23[make_pair(2,3)]+=1;
    //    dos23[make_pair(2,3)]+=1;
    //    dos23[make_pair(2,4)]+=1;
    //    dos23[make_pair(2,3)]+=1;
    //    dos23[make_pair(3,4)]+=1;


    //    for(auto& p1: dos23)
    //    {
    //        cout<<p1.first.first<< "\t";
    //        cout<<p1.first.second<< "\t";
    //        cout<<p1.second<<endl;
    //    }


    do {

        //cout<<sys2.state.toString()<<endl;
        e1=round(sys2.E()*10000);
        cl_v = check_closed_vertex(amount_of_circle,neighbours_for_cicle,is_vertex_closed,sys2,type_of_system);
        m_cl = init_calc_max_cluster(amount_of_circle,is_vertex_closed,neighbours_for_vertex);

        dos2_cl[make_pair(e1,cl_v)]+=2;
        dos2_mcl[make_pair(e1,m_cl)]+=2;



        if(ti%5368709==0){
            cout<<"Status: "<<ti/5368709<<endl;
        }

        ti++;
    } while (sys2.state.halfNext());
    //fq<<"!!!minimalE="<<minimalE;
    //fq<<"!!!maximumE="<<maximumE;

    for(auto& p1: dos2_cl)
    {
        fq<<p1.first.first/10000.<< "\t";
        fq<<p1.first.second<< "\t";
        fq<<p1.second<<endl;
    }
    cout<<"finish!"<<endl;

    for(auto& p1: dos2_mcl)
    {
        fqmcl<<p1.first.first/10000.<< "\t";
        fqmcl<<p1.first.second<< "\t";
        fqmcl<<p1.second<<endl;
    }
    cout<<"finish!"<<endl;
}

int main(int argc, char *argv[])
{
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;

    int CO = Prohod_MC_sampling/out_iterations;

    //srand(rank+time(NULL));

    //config::Instance()->srand(rank+time(NULL));
    config::Instance()->srand(1);

    //PartArray sys2;
    //sys2.load(examplename_mfsys);

    config::Instance()->set2D();
    SquareIsingLattice sys2;
    sys2.dropSquareLattice(4,4);
    sys2.setHamiltonian(hamiltonianIsing);


    //uniform_int_distribution<int> distr_int(0,sys2.size()-1);
    //uniform_real_distribution<double> distr_double(0.,1.);

    //    for(unsigned i=0;i<sys2.size();i++)
    //    {
    //        cout<<i<<"  "<<sys2[i]->pos.x<<"    "<<sys2[i]->pos.y<<endl;
    //    }

    double center_to_center_distance;   // расстояние между центрами
    double center_to_particle_distance;
    double particle_to_particle_nn_distance;

    if(type_of_system==0){              // для гексагональной
        center_to_center_distance=sqrt(3.);
        center_to_particle_distance=sqrt(3.)/2.;
        particle_to_particle_nn_distance=sqrt(3.)/2.;
    }
    else if(type_of_system==1){         // для треугольной
        center_to_center_distance=sqrt(3.)/3.;
        center_to_particle_distance=sqrt(3.)/6.;
        particle_to_particle_nn_distance=l/2.;
    }


    std::vector<double> center_of_circle;       // вектор хранящий центры гексагонов
    if(type_of_system==0){
        HC_create_center_of_circle(center_of_circle,sys2,l);    // ф-я, определяет параметр решетки l и вычисляет центры гексагонов, работает только с HC формой цветка,
    }
    else if (type_of_system==1) {
        Triangular_create_center_of_circle(center_of_circle,sys2,l);    // ф-я, определяет параметр решетки l и вычисляет центры треугольников
    }
    //    cout<<"Centers:"<<endl;
    //    for(unsigned i=0;i<center_of_circle.size();i+=2)
    //    {
    //        cout<<i<<"  "<<center_of_circle[i]<<"    "<<center_of_circle[i+1]<<endl;
    //    }

    int amount_of_circle = center_of_circle.size()/2;   // количество центров
    int amount_of_closed_vertex=0;  // количество замкнутых вертексов


    int neighbours_for_cicle[amount_of_circle][6];  // двумерный массив, хранящий 6 спинов каждого вертекса.
    int neighbours_for_vertex[amount_of_circle][6];  //двумерный  массив, хранящий макс. 6 соседних вертексов
    int is_vertex_closed[amount_of_circle];         // замкнут ли данный вертекс, 1-замкнут, 0-нет
    int Max_cluster_GS;                    // максимальный кластер

    //////////////////////// одноразовые команды///////////////////////////////

    if(type_of_system==0||type_of_system==1||type_of_system==2){
        ordering_vertex_spins(amount_of_circle,center_of_circle,neighbours_for_cicle,sys2,center_to_particle_distance,type_of_system); // заполнение neighbours_for_cicle и выстраивание очередности соседних спинов по кругу.

        for(int i=0;i<amount_of_circle;++i){
            for(int j=0;j<6;++j){
                neighbours_for_vertex[i][j]=-1;  // заполняем -1
            }
        }

        // определение соседних вертексов
        if(type_of_system==0){
            for(int i=0;i<amount_of_circle;i++){    // заполнение массива neighbours_for_vertex

                int jk=0;
                for(int j=0;j<amount_of_circle;j++){

                    if (g_isNeighbours_for_vertex(center_of_circle[2*i],center_of_circle[2*i+1],center_of_circle[2*j],center_of_circle[2*j+1],center_to_center_distance)){
                        neighbours_for_vertex[i][jk]=j;
                        jk++;
                    }
                    if(jk==6)
                        break;

                }
            }
        }
        else if(type_of_system==1){
            for(int i=0;i<amount_of_circle;i++){    // заполнение массива neighbours_for_vertex

                int jk=0;
                for(int j=0;j<amount_of_circle;j++){

                    if (g_isNeighbours_for_vertex(center_of_circle[2*i],center_of_circle[2*i+1],center_of_circle[2*j],center_of_circle[2*j+1],center_to_center_distance)){
                        neighbours_for_vertex[i][jk]=j;
                        jk++;
                    }
                    if(jk==3)
                        break;

                }
            }
        }


        //////////////////////////////////////////////////////////////////////

        //    cout<<"neighbours_for_vertex:"<<endl;
        //    for(int i=0;i<amount_of_circle;i++){
        //        cout<<i<<":  ";
        //        for(int j=0;j<6;j++){
        //            cout<<neighbours_for_vertex[i][j]<<"  ";
        //        }
        //        cout<<endl;
        //    }

        //    cout<<"neighbours_for_cicle:"<<endl;
        //    for(int i=0;i<amount_of_circle;i++){
        //        cout<<i<<":  ";
        //        for(int j=0;j<6;j++){
        //            cout<<neighbours_for_cicle[i][j]<<"   ";
        //        }
        //        cout<<endl;
        //    }


        amount_of_closed_vertex=check_closed_vertex(amount_of_circle,neighbours_for_cicle,is_vertex_closed,sys2,type_of_system); // пересчет is_vertex_closed // количество вихрей
        cout<<"amount_of_closed_vertex"<<amount_of_closed_vertex<<endl;
        // Внимание!! функция портит массив is_vertex_closed! Подсчет макс. кластера
        Max_cluster_GS=init_calc_max_cluster(amount_of_circle,is_vertex_closed,neighbours_for_vertex); // Внимание!! функция портит массив is_vertex_closed!
        cout<<"Max_cluster_GS"<<Max_cluster_GS<<endl;

        //    char *s1 = new char[sys2.size()];
        //    s1=strcpy(s1,sys2.state.toString().c_str());
        //    sys2.state.fromString(s1);
        //    delete[] s1;

        //exact_dos(is_vertex_closed, neighbours_for_cicle, neighbours_for_vertex, amount_of_circle, sys2);     //точное вычисление dos

        // Внимание!! функция портит массив is_vertex_closed! Подсчет макс. кластера

    }
    //////////// для вывода ///////////////

    /////////////////////////////////////////////////////////////////////////////////
    double *recvAE = new double[size];
    double *recvAE2 = new double[size];
    double *recvAE4 = new double[size];
    double *recvAM = new double[size];
    double *recvAPP1 = new double[size];
    double *recvAPP2 = new double[size];

    double *recvHeatCapacity = new double[size];
    double *recvVospr = new double[size];

    double *recvBCenergy = new double[size];
    double *recvBCmagn = new double[size];
    double *recvBCPP1 = new double[size];
    double *recvBCPP2 = new double[size];


    ofstream outAPP1("APP1.dat",ios::out);            //средний параметр порядка
    ofstream outAPP2("APP2.dat",ios::out);          //средний параметр порядка2
    ofstream outAE("energyAVG.dat",ios::out);       //средняя энергия
    ofstream outAE2("energyAVG2.dat",ios::out);     //средняя энергия квадрат
    ofstream outAE4("energyAVG4.dat",ios::out);     //средняя энергия 4 степень
    ofstream outAM("AMagn.dat",ios::out);           //Средняя намагниченность
    ofstream outHeatCapacity("C.dat",ios::out);     //Теплоемкость
    ofstream outVospr("X.dat",ios::out);            //Восприимчивость


    // биндеры
    ofstream outBCenergy("BC_energy.dat",ios::out); //биндер по энергии
    ofstream outBCmagn("BC_magn.dat",ios::out);     //биндер по намагниченности
    ofstream outBCPP1("BC_PP1.dat",ios::out);         //биндер по параметру порядка
    ofstream outBCPP2("BC_PP2.dat",ios::out);       //биндер по параметру порядка2
    /////////////////////////////////////////////////////////////////////////////////////

    ////////// для вывода всех текущих значений до усреднения///////

    //if(out_current_values){
    double *recvE;
    double *recvPP1;
    double *recvPP2;
    ofstream outE("energy.dat",ios::out);
    ofstream outPP1("PP1.dat",ios::out);
    ofstream outPP2("PP2.dat",ios::out);

    vector<double> CValue;          // E current value
    vector<double> CValuePP1;       // PP current value
    vector<double> CValuePP2;       // PPF current value


        recvE = new double[out_iterations*size];
        recvPP1 = new double[out_iterations*size];
        recvPP2 = new double[out_iterations*size];

    //////////////////////////////////


    /////////////////////////MC/////////////////////////

    double temperature=0;
    double log_maxtemp=log(maxtemp);
    double log_mintemp=log(mintemp);
    double d_log_t = (log_maxtemp-log_mintemp)/size;

    temperature=set_temp(rank,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t);


    double RBCenergy; // результирующий БК энергии
    double RBCmagn; // результирующий БК магнитный
    double RBCPP2; // результирующий БК параметра порядка 2
    double RBCPP1; // результирующий БК параметра порядка

    double SumenergyVar,SumenergyVar2, SumenergyVar4;
    double temMagn,MagnVar,MagnVar2, MagnVar4;
    double Aenergy=0;
    double Aenergy2=0;
    double Aenegry4=0;
    double Amagn2=0;
    double Amagn4=0;
    double X_o=0;



    double e1=sys2.E();
    double pp1=0;
    double pp2=0;
    double pp1_2=0;
    double pp2_2=0;
    double pp1_4=0;
    double pp2_4=0;
    double app1=0;
    double app2=0;
    double app1_2=0;
    double app2_2=0;
    double app1_4=0;
    double app2_4=0;


    // для вывода C

    double C_o=0;
    //double *recvHeatCapacity = new double[size];
    //ofstream outHeatCapacity("C.dat",ios::out);

    // для репличного обмена
    double exchange_t=0;
    double exchange_e=0;
    char exchange_state[sys2.size()];
    int yes_no_exchange=0;
    double probability_of_exchange=0;

    unsigned long long Prohod=0;
    unsigned long long o_index=0;

    ////////////Monte-Carlo разогрев////////////////////
    //MC_step(sys2,Prohod_MC_equilibration,e1,temperature);
    //////////// MC steps //////////////////////////////
    while(Prohod<Prohod_MC_sampling){
        e1=sys2.E();
        //обмен конфигурациями каждые between_exchange шагов///////////////////////////////
        if(replica_on){
            if(Prohod%between_exchange==0){
                for(int ii=size-1;ii>0;--ii){
                    if(rank==ii){
                        MPI_Send(&temperature,1, MPI_DOUBLE, ii-1, 0, MPI_COMM_WORLD);
                        MPI_Send(&e1,1, MPI_DOUBLE, ii-1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&yes_no_exchange,1, MPI_INT, ii-1, 0, MPI_COMM_WORLD,&status);
                        if(yes_no_exchange==0){}

                        else{
                            MPI_Send(sys2.state.toString().c_str(),sys2.size(), MPI_CHAR, ii-1, 0, MPI_COMM_WORLD);
                            MPI_Recv(&exchange_state,sys2.size(), MPI_CHAR, ii-1, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&e1,1, MPI_DOUBLE, ii-1, 0, MPI_COMM_WORLD,&status);
                            sys2.state.fromString(exchange_state);
                            //e1=sys2.E();
                        }

                    }
                    if(rank==ii-1){
                        MPI_Recv(&exchange_t,1, MPI_DOUBLE, ii, 0, MPI_COMM_WORLD,&status);
                        MPI_Recv(&exchange_e,1, MPI_DOUBLE, ii, 0,MPI_COMM_WORLD,&status);
                        probability_of_exchange=exp(((1/temperature)-(1/exchange_t))*(e1-exchange_e));
                        if (probability_of_exchange<rand()/RAND_MAX){
                            yes_no_exchange=0;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, ii, 0, MPI_COMM_WORLD);

                        }
                        else{
                            yes_no_exchange=1;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, ii, 0, MPI_COMM_WORLD);

                            MPI_Send(sys2.state.toString().c_str(),sys2.size(), MPI_CHAR, ii, 0, MPI_COMM_WORLD);
                            MPI_Recv(&exchange_state,sys2.size(), MPI_CHAR, ii, 0, MPI_COMM_WORLD,&status);
                            MPI_Send(&e1,1, MPI_DOUBLE, ii, 0, MPI_COMM_WORLD);
                            e1=exchange_e;
                            sys2.state.fromString(exchange_state);
                            //e1=sys2.E();


                        }

                    }

                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        //////////////////iterations_wo_out//////////////////////
        e1=sys2.E();
        MC_step(sys2,iterations_wo_out,e1,temperature);

        ////////////////////out////////////////////
        if(Prohod%CO==0){
            o_index++;

            //e1=sys2.E();

            SumenergyVar = e1;     /////////////
            SumenergyVar2=SumenergyVar*SumenergyVar;
            SumenergyVar4=SumenergyVar*SumenergyVar*SumenergyVar*SumenergyVar;
            Aenergy=(SumenergyVar+(o_index-1)*Aenergy)/(o_index);
            Aenergy2=(SumenergyVar2+(o_index-1)*Aenergy2)/(o_index);
            Aenegry4=(SumenergyVar4+(o_index-1)*Aenegry4)/(o_index);



            MagnVar=sys2.M().x;		/////////////
            MagnVar2=MagnVar*MagnVar;
            MagnVar4=MagnVar*MagnVar*MagnVar*MagnVar;
            temMagn=(MagnVar+(o_index-1)*temMagn)/(o_index);   // !!!*N куммулятив смешанного кластера
            Amagn2=(MagnVar2+(o_index-1)*Amagn2)/(o_index);
            Amagn4=(MagnVar4+(o_index-1)*Amagn4)/(o_index);

            pp1=(double) check_closed_vertex(amount_of_circle,neighbours_for_cicle,is_vertex_closed,sys2,type_of_system)/(double)Max_cluster_GS;
            pp1_2=pp1*pp1;
            pp1_4=pp1*pp1*pp1*pp1;
            app1=(pp1+(o_index-1)*app1)/(o_index);
            app1_2=(pp1_2+(o_index-1)*app1_2)/(o_index);
            app1_4=(pp1_4+(o_index-1)*app1_4)/(o_index);


            pp2=(double)init_calc_max_cluster(amount_of_circle,is_vertex_closed,neighbours_for_vertex)/(double)Max_cluster_GS;
            pp2_2=pp2*pp2;
            pp2_4=pp2*pp2*pp2*pp2;
            app2=(pp2+(o_index-1)*app2)/(o_index);
            app2_2=(pp2_2+(o_index-1)*app2_2)/(o_index);
            app2_4=(pp2_4+(o_index-1)*app2_4)/(o_index);



                if(out_current_values){
                    CValue.insert(CValue.end(),e1);
                    CValuePP1.insert(CValuePP1.end(),pp1);
                    CValuePP2.insert(CValuePP2.end(),pp2);
                }



        }



        if(rank==0)
        {
            if(Prohod%(Prohod_MC_sampling/100)==0)
                cout<<"Status: "<<Prohod/(Prohod_MC_sampling/100)<<endl;
        }
        Prohod++;

    }

    C_o=((Aenergy2)-(Aenergy*Aenergy))/(temperature*temperature*sys2.size()*sys2.size()); //  теплоёмкость //  теплоёмкость
    X_o=((Amagn2)-(temMagn*temMagn))/(temperature*temperature*sys2.size()*sys2.size());
    RBCenergy = 1 -(Aenegry4/(3*pow(Aenergy2,2)));  // результирующий БК энергии
    RBCmagn = 1-(Amagn4/(3*pow(Amagn2,2)));  // результирующий БК магнитный
    RBCPP1 = 1-(app1_4/(3*pow(app1_2,2)));  // результирующий БК параметра порядка 2
    RBCPP2 = 1-(app2_4/(3*pow(app2_2,2)));  // результирующий БК параметра порядка


    if(out_current_values)
    {
        MPI_Gather(&CValue[0], out_iterations, MPI_DOUBLE, recvE, out_iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&CValuePP1[0], out_iterations, MPI_DOUBLE, recvPP1, out_iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&CValuePP2[0], out_iterations, MPI_DOUBLE, recvPP2, out_iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }

    //MPI_Gather(&e_AVG, 1, MPI_DOUBLE, recvAE, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(&PP1_AVG, 1, MPI_DOUBLE, recvAPP1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(&PP2_AVG, 1, MPI_DOUBLE, recvAPP2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(&RBC, 1, MPI_DOUBLE, recvBC, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(&RBC_PP1, 1, MPI_DOUBLE, recvBC_PP1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(&C_o, 1, MPI_DOUBLE, recvHeatCapacity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gather(&X_o, 1, MPI_DOUBLE, recvVospr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&C_o, 1, MPI_DOUBLE, recvHeatCapacity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Теплоемкость
    MPI_Gather(&RBCenergy, 1, MPI_DOUBLE, recvBCenergy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК энергию
    MPI_Gather(&RBCmagn, 1, MPI_DOUBLE, recvBCmagn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК магн
    MPI_Gather(&RBCPP1, 1, MPI_DOUBLE, recvBCPP1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PP1
    MPI_Gather(&RBCPP2, 1, MPI_DOUBLE, recvBCPP2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PP2
    MPI_Gather(&Aenergy, 1, MPI_DOUBLE, recvAE, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenergy2, 1, MPI_DOUBLE, recvAE2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenegry4, 1, MPI_DOUBLE, recvAE4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temMagn, 1, MPI_DOUBLE, recvAM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&pp1, 1, MPI_DOUBLE, recvAPP1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&pp2, 1, MPI_DOUBLE, recvAPP2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if(rank==0)
    {
        //            // out info
        //            outAE<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            outAPP1<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            outAPP2<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            outBCPP2<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            outBCPP1<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            outHeatCapacity<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            if(out_current_values){
        //                outE<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //                outPP1<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //                outPP2<<"#Info:"<<endl<<"#Numers of particles: "<<sys2.size()<<endl<<"#Numbers of circle: "<<amount_of_circle<<endl<<"#Number of vorteces:"<<amount_of_closed_vertex<<endl<<"#Max_cluster: "<<Max_cluster_GS<<endl<<"#Min T: "<<mintemp<<",   Max T: "<<maxtemp<<endl<<"#Replica exchage: "<<replica_on<<",  Out current values:"<<out_current_values<<endl<<"#First heating (Numbers of MC steps before start):"<<first_iteration_wo_count<<endl<<"#Number of output averaged MC steps: "<<out_iterations<<endl<<"#Number of MC steps between replica exchange"<<between_exchange<<endl<<"#Number of replica exchanges"<<out_iterations/between_exchange<<endl<<"#Number of iterations without output(Heating if each MC out steps): "<<iterations_wo_out<<endl<<"#Total number of MC steps"<<TOTAL_MC_Steps<<endl;
        //            }

        for(int i=0; i<size; ++i)
        {
            if(out_current_values){
                for(unsigned long long j=0; j<out_iterations; ++j)
                {
                    outE<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvE[i*out_iterations+j]<<endl;
                    //outPP1<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvPP1[i*out_iterations+j]<<endl;
                    //outPP2<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvPP2[i*out_iterations+j]<<endl;
                }
            }
            //outAE<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvAE[i]/(out_iterations)<<endl;
            //                outAPP1<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvAPP1[i]/(out_iterations)<<endl;
            //                outAPP2<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvAPP2[i]/(out_iterations)<<endl;
            //                outBCPP2<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvBC[i]<<endl;
            //                outBCPP1<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvBC_PP1[i]<<endl;
            //outHeatCapacity<<set_temp(i,temperature_distribution,maxtemp, mintemp,size,log_mintemp,d_log_t)<<'\t'<<recvHeatCapacity[i]<<endl;
            outAPP1<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPP1[i]<<endl;
            outAPP2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPP2[i]<<endl;
            outAM<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM[i]<<endl;
            outAE<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE[i]<<endl;
            outAE2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE2[i]<<endl;
            outAE4<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE4[i]<<endl;

            outHeatCapacity<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvHeatCapacity[i]<<endl;                 // Теплоемкость

            outVospr<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvVospr[i]<<endl;  //воспириимчивость
            outBCenergy<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCenergy[i]<<endl; // БК energy
            outBCmagn<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCmagn[i]<<endl; // БК magn
            outBCPP1<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPP1[i]<<endl; // БК PPF
            outBCPP2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPP2[i]<<endl; // БК PP
        }
    }


    //if(out_current_values)
        delete [] recvE;
        delete [] recvPP1;
        delete [] recvPP2;


    delete []recvAPP1;
    delete []recvAPP2;
    delete []recvAM;
    delete []recvAE;
    delete []recvAE2;
    delete []recvAE4;
    delete []recvBCenergy;
    delete []recvBCmagn;
    delete []recvBCPP1;
    delete []recvBCPP2;
    delete []recvHeatCapacity;
    delete []recvVospr;

    MPI_Finalize();

    cout << "finish\n";
    return 0;

}
