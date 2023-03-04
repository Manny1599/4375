/*
Logistic Regression for Titanic Data
Manuel Romero
Sagar Darnal

*/
using namespace std;
#include <iostream>;
#include <fstream>;
#include <vector>;
#include <algorithm>;
#include <math.h>;
#include <cmath>
#include <chrono>


double e = 2.71828;

double sigmoid(double z)
{
    return (1 / (1 + exp(-z)));
}

//a function to find the mean of a numeric vector
//takes in a vector returns a double, size of vector, what int to look for, survived vector
double vecMean(vector<double> vect, int n, int b, vector<double> vec)
{
    double mean = 0, tot = 0;
    for(int i =0; i < n; i++)//train data mean
    {
        if(vec[i] == b)//if survived or died given b
        {
            mean += vect[i];
            tot++;//total number of survuved or dead
        }

    }
    mean /= tot;
    return mean;
}
//-----------------------------------------------------------------------


//a function that computes the variance of a vector
//takes in a vector, size of vecotor and what int to find, survived vector
double vecVar(vector<double> vect, int n, int b, vector<double> vec)
{
   double var = 0.0, vectmean, tot=0;
   vectmean = vecMean(vect, n,b,vec);

   for (int i = 0; i < n; i++)
   {
       if(vec[i] == b)//survived or dead given b
       {
           tot++;//getting total of survivors or dead given b
       }
    }
   for (int i = 0; i < n; i++)
   {
       if(vec[i] == b)//survived or dead given b
       {
           var += (pow((vect[i] - vectmean), 2)) / (tot-1);
       }
    }
   return var;

}//end vecVar
//-----------------------------------------------------------------------

//function to display accuracy, sensitivity and specificity of a logistic function
//takes in intercept, slope, what number the test data starts at, predictor vector and the actual values vector(test data), and size of the data(end of test data)
//returns nothing
void accuracy(double b, double m, int test, vector<double> x, vector<double> y, int sizex)
{
    double tp = 0;
    double fp = 0;
    double tn = 0;
    double fn = 0;

    double predicted;
    for(int i = test; i < sizex; i++)
    {

        predicted = b + m * x[test];
        if(predicted < .5)//making predictions 0 or 1
        {
            predicted = 0.0;//died
        }
        else
        {
            predicted = 1.0;//survived
        }


        if(predicted == y[test] && predicted == 1.0) // tp
        {
            tp++;
        }
        else if (predicted != y[test] && predicted == 1.0)//fp
        {
            fp++;
        }
        else if (predicted == y[test] && predicted == 0.0)//tn
        {
            tn++;
        }
        else if (predicted != y[test] && predicted == 0.0)//fn
        {
            fn++;
        }
        else
        {
            cout << "Erorr\n";
        }
        //cout<< predicted << "    " << y[test]<< "\n" ;
        test ++;
    }
    cout << "tp: " << tp;
    cout << "\nfp: " << fp;
    cout << "\ntn: " << tn;
    cout << "\nfn: " << fn;

    cout << "\nAccuracy is  :" <<  (tp + tn)/(tp+fp+fn+tn) << "\n";
    cout << "Sensitivity is  :" <<  tp/(tp+fn) << "\n";
    cout << "Specificity is  :" <<  tn/(tn+fp) << "\n";


}

//calulating age likelihood given survival or death
double calcAge(double var, double mean, double age)
{
    return ((1/ (2*M_PI * var))*exp(-(pow((age - mean),2)/(2*var))));
}


//-----------------------------------------------------------------------
int main(int argc, char** argv)
{
    const int MAX_LEN = 10000;
    vector<double> sexc(MAX_LEN); //sex
    vector<double> pc(MAX_LEN); //pclass vector
    vector<double> surv(MAX_LEN); //survive vector
    vector<double> age_p(MAX_LEN); //age vecrtor
    //code given by Proffesor on reading from csv file
    ifstream inFS; //input file stream
    string line;
    string col1, pclass, survived, sex, age;


     cout << "Opening file titanic_project.csv" << endl;
    inFS.open("/Users/fernandoromero/Documents/titanic_project.csv");
    if(!inFS.is_open())
    {
        cout << "Could not open file titanic_project.csv" << endl;
        return 1;
    }

    //file is open and we can work on it
    //titanic_procect.csv has 5 int columns

    cout << "Reading Line 1" << endl;
    getline(inFS, line);
    //displaying heading
    cout << "HEADING: " << line << endl;
    //timing training time
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    int numObservations = 0;
    while (inFS.good())
    {
        getline(inFS, col1, ',');
        getline(inFS, pclass, ',');
        getline(inFS, survived, ',');
        getline(inFS, sex, ',');
        getline(inFS, age, '\n');


        pc.at(numObservations) = stof(pclass);
        surv.at(numObservations) = stof(survived);
        sexc.at(numObservations) = stof(sex);
        age_p.at(numObservations) = stof(age);

        numObservations++;
    }

    pc.resize(numObservations);
    surv.resize(numObservations);
    sexc.resize(numObservations);
    age_p.resize(numObservations);

    cout << "New Length : " << pc.size() <<endl;

    cout << "Closing file titanic_project.csv" << endl;
    inFS.close(); //Done with the file

    // ------------------------------------------------------------------ data read ---------------------------------------
    //doing 500,000 itereations
    auto t1 = high_resolution_clock::now(); //start time

    double weight[] = {0.0, 0.0,0.0,0.0};//class , survived, sex ,age
    double learningRate = .001;

    for(int j =0; j<500000; j++)
    {
        double error1=0, error2=0;//error 1 for sex and error 2 for actual

        for(int i=0; i<800; i++)//data train
        {
            double p = sexc[i] * weight[2] + weight[1];
            double h = sigmoid(p);
            error1 += (h - surv[i])*sexc[i] ;
            error2 += h - surv[i];
        }
        weight[1] = weight[1] - learningRate * error2/800;//intercept
        weight[2] = weight[2] - learningRate * error1/800;//slope

    }

    cout << "Intercept :  " << weight[1] << "   slope :   " << weight[2] << "\n";
    auto t2 = high_resolution_clock::now(); //end time

    accuracy(weight[1], weight[2], 800, sexc, surv, pc.size());
    duration<double, std::milli> ms_double = t2 - t1;
    cout <<"Training Time of algorithm Logistic Regerssion: " <<  ms_double.count() << "ms\n";

    //
    // ------------------------------------------------------------------ Naive Bayes -------------------------------------
    //prior prob
    auto t12 = high_resolution_clock::now(); //start time
    double survival = 0, died = 0, tot = 0;
    int survn = 0, diedn = 0;
    for(int i = 0; i < 800; i++) // train data
    {
        if(surv[i] == 1)
        {
            survival++;
            survn++;
        }
        else
        {
            died++;
            diedn++;
        }
        tot++;

    }

    died = died/tot; //prior prob of died
    survival = survival/tot; //prior prob of surv
    cout << "\nThe prior probability of survived = no(0) : " << died;
    cout << "\nThe prior probability of survived = yes(1) : " << survival;


    //conditional prob
    //lilihood for classes
    //classA survived given class 1, classAT number of memebers in class 1 , same for class 2 and 3
    double class1 = 0, class2 = 0, class3 = 0, class1T = 0, class2T = 0, class3T = 0;
    double l1, l2, l3, l1D, l2D, l3D;//likelihood survived and died given class
    for(int i = 0; i < 800; i++) // train data
    {
        if( pc[i] == 1)// if in class 1
        {
            class1T++;
            if(surv[i] == 1) // survived
            {
                class1++;
            }
        }//end of class 1

        else if( pc[i] == 2)// if in class 2
        {
            class2T++;
            if(surv[i] == 1) // survived
            {
                class2++;
            }
        }//end of class 2

        else //( pc[i] == 3)// if in class 3
        {
            class3T++;
            if(surv[i] == 1) // survived
            {
                class3++;
            }
        }//end of class 3
    }//end of data seperation
    //survived likelihood for classes

    l1 = class1/survn; //class 1 survival divided by tot survival
    l2 = class2/survn; // class 2 survival divided by tot survival
    l3 = class3/survn; // class 3 survival divided by tot survival

    //died likelihood for classes
    l1D = (class1T - class1) /diedn; //class 1 died divided by tot died
    l2D = (class2T - class2) /diedn; //class 2 died divided by tot died
    l3D = (class3T - class3) /diedn; //class 3 died divided by tot died

    cout << "\nLikelihoods for survival given class: \nclass \tsurvived \tdied\nc1 : \t" <<l1 <<"\t" << l1D;
    cout << "\nc2 : \t" <<l2 <<"\t" << l2D;
    cout << "\nc3 : \t" <<l3 <<"\t" << l3D << "\n";

    //died likelihood for sex assuming 0 is for female
    double fsurv = 0, msurv = 0, ftot = 0, mtot = 0, lf, lm, lfd, lmd;
    for(int i = 0; i < 800; i++) // train data
    {
        if( sexc[i] == 0)// if female/woman
        {
            ftot++;
            if(surv[i] == 1) // survived
            {
                fsurv++;
            }
        }//end of female

        else //( sexc[i] == 1)// if male
        {
            mtot++;
            if(surv[i] == 1) // survived
            {
                msurv++;
            }
        }//end of male
    }//end of data seperation

    //calculating likelihood for male and female
    lf = fsurv/survn; //female survival divided by tot survival
    lm = msurv/survn; // male survival divided by tot survival

    //died likelihood for classes
    lfd = (ftot - fsurv) /diedn; //female died divided by tot died
    lmd = (mtot - msurv) /diedn; //male died divided by tot died

    //print results
    cout << "Likelihoods for survival given class: \nSex \t\tsurvived \tdied\nfemale : \t" <<lf <<"\t" << lfd;
    cout << "\nmale : \t\t" <<lm <<"\t" << lmd << "\n";

    //end of likelihood given sex


    //calculating likelihood given age
    //getting mean and variance of age
    double ageM, ageV, ageMD, ageVD;
    ageM = vecMean(age_p, 800, 1, surv);//mean for train data for survived
    ageV = vecVar(age_p, 800, 1, surv);//variance for survived
    ageMD = vecMean(age_p, 800, 0, surv);//mean for train data for dead
    ageVD = vecVar(age_p, 800, 0 , surv);// variance for dead

    cout << "\nsurvived\tage mean\tage varr\n";
    cout << "T\t\t" << ageM <<"\t\t" << ageV;
    cout << "\nF\t\t" << ageMD <<"\t\t" << ageVD << "\n";

    auto t22 = high_resolution_clock::now(); //end time


    duration<double, std::milli> ms_double2 = t22 - t12;
    cout <<"Training Time of algorithm Naive Bayes: " <<  ms_double2.count() << "ms\n";



    double num_s, num_d,tp=0,tn=0,fp=0,fn=0, predicted, predictedD, denom;

    //runing test
    for(int i = 800; i < 1024; i++)
    {
        //class * sex * age * prior
        if(pc[i] ==1)//if in first class
        {
            if(sexc[i] ==0 )//if female
            {
                //likelihood of survival or death given they were in first class, female, and age x
                predicted = l1 * lf * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l1D * lfd * calcAge(ageVD,ageMD,age_p[i]) * died;
            }
            else//if male
            {
                //likelihood of survival given they were in first class, male, and age x
                predicted = l1 * lm * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l1D * lmd * calcAge(ageVD,ageMD,age_p[i]) * died;

            }
        }
        if(pc[i] ==2)//if in second class
        {
            if(sexc[i] ==0 )//if female
            {
                //likelihood of survival or death given they were in 2nd class, female, and age x
                predicted = l2 * lf * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l2D * lfd * calcAge(ageVD,ageMD,age_p[i]) * died;
            }
            else//if male
            {
                //likelihood of survival given they were in 2nd class, male, and age x
                predicted = l2 * lm * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l2D * lmd * calcAge(ageVD,ageMD,age_p[i]) * died;

            }
        }
        if(pc[i] ==3)//if in 3rd class
        {
            if(sexc[i] ==0 )//if female
            {
                //likelihood of survival or death given they were in 3rd class, female, and age x
                predicted = l3 * lf * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l3D * lfd * calcAge(ageVD,ageMD,age_p[i]) * died;
            }
            else//if male
            {
                //likelihood of survival given they were in 3rd class, male, and age x
                predicted = l3 * lm * calcAge(ageV,ageM,age_p[i]) * survival;
                predictedD = l3D * lmd * calcAge(ageVD,ageMD,age_p[i]) * died;

            }
        }


        denom = predicted + predictedD;
        predicted = predicted/denom;

        //got predicted values now
        if(predicted < .5) //rounding and comparing
        {
            predicted = 0;
        }
        else
        {
            predicted = 1;
        }

         //we will use only predicted since that is the prediction for survival
         if(predicted == surv[i] && predicted == 1.0) // tp
        {
            tp++;
        }
        else if (predicted != surv[i] && predicted == 1.0)//fp
        {
            fp++;
        }
        else if (predicted == surv[i] && predicted == 0.0)//tn
        {
            tn++;
        }
        else if (predicted != surv[i] && predicted == 0.0)//fn
        {
            fn++;
        }
        else
        {
            cout << "Erorr\n";
        }
    }
        //cout<< predicted << "    " << y[test]<< "\n" ;
    cout << "tp: " << tp;
    cout << "\nfp: " << fp;
    cout << "\ntn: " << tn;
    cout << "\nfn: " << fn;

    cout << "\nAccuracy is  :" <<  (tp + tn)/(tp+fp+fn+tn) << "\n";
    cout << "Sensitivity is  :" <<  tp/(tp+fn) << "\n";
    cout << "Specificity is  :" <<  tn/(tn+fp) << "\n";


    return 0;
}















