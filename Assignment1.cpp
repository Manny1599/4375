//Data Exploration
//Manuel Romero
//CS 4375.004
using namespace std;
#include <iostream>;
#include <fstream>;
#include <vector>;
#include <algorithm>;
#include <math.h>;

//a function to find the sum of a numeric vector
//takes in a vector returns a double
double vecSum(vector<double> vect)
{
    double sum = 0;
    for (int i = 0; i < vect.size(); i++)
        sum +=vect[i];
    return sum;

}//end vecSum
//-----------------------------------------------------------------------


//a function to find the mean of a numeric vector
//takes in a vector returns a double
double vecMean(vector<double> vect)
{
    return (vecSum(vect)/ vect.size());
}//end vecMean
//-----------------------------------------------------------------------

//a function to find the median of a numeric vector
//takes in a vector returns a double
double vecMedian(vector<double> vect)
{
    int middle = vect.size()/2;
    double median;
    //sorting the vector in ascending order
    sort(vect.begin(), vect.end());
    //if even getting the average of the two numbers in the middle
    //else just getting the middle number
    if (vect.size() % 2 == 0)
    {
        median = (vect[middle] + vect[(middle+1)])/2;
        return median;
    }
    else
    {
        return median = vect[middle];
    }
}//end vecMedian
//-----------------------------------------------------------------------

//a function to find the range of a numeric vector
//takes in a vector returns a double
double vecRange(vector<double> vect)
{
    double range;
    //sorting the vector in ascending order
    sort(vect.begin(), vect.end());
    return range = (vect[(vect.size()-1)] - vect[0]);
}//end vecRange
//-----------------------------------------------------------------------

//a function to compute covariance between rm and medv
//takes in two vectors returns a double
double vecCov(vector<double> vectx, vector<double> vecty)
{
    double Cov = 0, xmean, ymean;
    int n = vectx.size();//size of vector

    //getting the means
    xmean = vecMean(vectx);
    ymean = vecMean(vecty);

    //calculating the covariance
    for (int i = 0; i < n; i++)
        Cov += ((vectx[i] - xmean) * ( vecty[i] - ymean)) / (n-1);
    return Cov;
}//end vecCov
//-----------------------------------------------------------------------


//a function that computes the sigma of a vector
//takes in a vector
double vecSigma(vector<double> vect)
{
   double sigma = 0.0, vectmean;
   double n = vect.size();//size of vector
   vectmean = vecMean(vect);
   for (int i = 0; i < n; i++)
   {
    sigma += (pow((vect[i] - vectmean), 2)) / (n-1);
   }
   sigma = sqrt(sigma);
   return sigma;

}//end vecSigma
//-----------------------------------------------------------------------

//a function to compute correlation between rm and medv
//takes in two vectors returns a double
double vecCorr(vector<double> vectx, vector<double> vecty)
{
    double corr, sigmaX, sigmaY;
    sigmaX = vecSigma(vectx);
    sigmaY = vecSigma(vecty);

    //sorting the vector in ascending order
    corr = (vecCov(vectx, vecty))/ (sigmaX * sigmaY);
    return corr;
}//end vecCorr
//-----------------------------------------------------------------------


//a function that prints out the stats for given vector
//takes in a vector
//prints sum, mean, median, range and covariance and correlation
void print_stats(vector<double> vect)
{
   cout << "The sum of the vector data is :      " << vecSum(vect) << endl;
   cout << "The mean of the vector data is :     " << vecMean(vect) << endl;
   cout << "The median of the vector data is :   " << vecMedian(vect) << endl;
   cout << "The range of the vector data is :    " << vecRange(vect) << endl;
}//end vecSigma
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
int main(int argc, char** argv)
{
    //code given by Proffesor on reading from csv file
    ifstream inFS; //input file stream
    string line;
    string rm_in, medv_in;
    const int MAX_LEN = 1000;
    vector<double> rm(MAX_LEN);
    vector<double> medv(MAX_LEN);

    //trying to open file
    cout << "Opening file Boston.csv" << endl;
    inFS.open("/Users/fernandoromero/Documents/HW1/Boston.csv");
    //inFS.open("Boston.csv");
    if(!inFS.is_open())
    {
        cout << "Could not open file Boston.csv" << endl;
        return 1;
    }
    //file is open and we can work on it
    //Boston.csv has two double columns
    cout << "Reading Line 1" << endl;
    getline(inFS, line);
    //displaying heading
    cout << "HEADING: " << line << endl;

    int numObservations = 0;
    while (inFS.good())
    {
        getline(inFS, rm_in, ',');
        getline(inFS, medv_in, '\n');

        rm.at(numObservations) = stof(rm_in);
        medv.at(numObservations) = stof(medv_in);
        numObservations++;
    }

    rm.resize(numObservations);
    medv.resize(numObservations);
    cout << "New Length : " << rm.size() <<endl;

    cout << "Closing file Boston.csv" << endl;
    inFS.close(); //Done with the file

    cout<< "Number of records: " << numObservations << endl;
    cout << "\nStats for rm" << endl;
    print_stats(rm);

    cout << "\nStats for medv" << endl;
    print_stats(medv);

    cout << "\n Covariance = " << vecCov(rm, medv) << endl;
    cout << "\n Corralation = " << vecCorr(rm, medv) << endl;
    cout << "\nProgram terminated.";

    return 0;
}//end main
