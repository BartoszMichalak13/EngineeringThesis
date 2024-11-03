// https://bitbucket.org/NachiketDan/dreyfus-wagner/src/master/

#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <cmath>
#include <queue>
#include <vector>
#include <memory>

using namespace std;

#define max 10

int p[max][max];
/*
 * All Pairs Shortest Path using Floyd's Algorithm
 */
void allpairshort(int a[max][max], int n)
{
  int k, i, j;
  for (k = 0; k < n; k++)
  {
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        if (a[i][k] + a[k][j] < a[i][j])
        {
          a[i][j] = a[i][k] + a[k][j];
          p[i][j] = k;
        }
      }
    }
  }
}

/*
 * Storing the shortest path
 */
void shortest(int i, int j)
{
  int k = p[i][j];
  if (k > 0)
  {
    shortest(i, k);
    cout << "  " << k << "  ";
    shortest(k, j);
  }
}

/*hashTable
Key = {1,2}  Value = 0
Key = {2,3}  Value = 1
{1,4}
{1,2,3}
....
...
{1,2,3,4}
{1,2,3,4} = {1,2}, {2,3}, {3,4}, {1,4}, {2,4}, {1,3} {1,2,3,4}
Target subset size = 3
E = 4
{1,2,3}
{2,3,4}
numberOfSubsets of size r = n!/(n-r)!*r!*/

map< set <int>, int> subsetMap;

vector< set< int > > generateSubsets(vector<int> C, int n,int graphsize)
{
  int count = pow(2, n);

  vector<set <int> > allSubsets;
  int subsetIndex = graphsize;
  // The outer for loop will run 2^n times to print all subset .
  // Here variable i will act as a binary counter

  for (int i = 0; i < count; i++)
  {
    // The inner for loop will run n times , As the maximum number of elements a set can have is n
    // This loop will generate a subset
    set<int> newSubset;
    for (int j = 0; j < n; j++)
    {
      // This if condition will check if jth bit in binary representation of  i  is set or not
      // if the value of (i & (1 << j)) is greater than 0 , include arr[j] in the current subset
      // otherwise exclude arr[j]
      if ((i & (1 << j)) > 0)
      {
          newSubset.insert(C[j]);
      }
    }

    if (newSubset.size() > 1)
    {
      allSubsets.push_back(newSubset);
      //subsetMap[newSubset] = subsetIndex;
      subsetMap.insert(make_pair(newSubset,subsetIndex));
      cout << "Adding subset with index" << subsetIndex << "\n";
      subsetIndex++;
    }
  }

  for (auto const &x : subsetMap)
  {
    cout << "Subset = ";
    set<int> element = x.first;
    for (auto f : element)
    {
      cout << f << " ";
    }

    cout << "\n";
    std::cout << "Index = " << x.second << std::endl;
  }

  cout<<"Completed subsets";

  return allSubsets;
}

void initializeSteinerTable(vector<vector<int>> &steinerDistanceTable, vector<vector<int>> &graph, int cSize)
{
  for (int i = 0; i < graph.size(); i++)
  {
    for (int j = 0; j < graph.size(); j++)
    {
      steinerDistanceTable[i][j] = graph[i][j];
    }
  }
}

int fetchIndexofMapofSets(set<int> subset)
{
  int key;
  for(map< set<int> , int >::const_iterator it = subsetMap.begin(); it != subsetMap.end(); ++it)
  {
    if (it -> first == subset)
    {
      key = it ->second;
    }
  }

  return key;
}

int calculateSteiner(vector<int> C, vector<vector<int>> graph, vector<set<int>> allSubsets,int q)
{
  // Remember that C = Y - {any one element} where Y = {Steiner terminals}
  // Key is index of D and value is distance
  //nodes = {1,2,3,4,5,6,7}
  //D = { {2,3}, {2,4}, {3,4}, {2,3,4} }
  int sizetable = allSubsets.size() + graph.size();
  vector< vector<int> > steinerDistanceTable(sizetable, vector<int> (graph.size()));

  //initialise steiner table
  initializeSteinerTable(steinerDistanceTable, graph, sizetable);

  int indexOfSubset = 0;

  set<set<int>>::iterator setIterator; // iterator for the "outer" structure

  for (int m = 2; m < C.size(); m++)
  {
    // for each subset
    for (int i = 0; i < allSubsets.size(); i++)
    {
      set<int> elementOfD = allSubsets[i];
      if(elementOfD.size() == m)
      {
        indexOfSubset = subsetMap[elementOfD];
        // initialise distance from elementOfD to every node in graph
        for (int i = 0; i < graph.size(); i++)
        {
          steinerDistanceTable[indexOfSubset][i] = 999;
        }

        // for every element in graph
        for (int j = 0; j < graph.size(); j++)
        {
          int u = 999 ;

          set<int>::iterator nodeIterator;
          for (nodeIterator = elementOfD.begin(); nodeIterator != elementOfD.end(); nodeIterator++)
          {
            // E is being considered as a single number from the set named elementOfD
            int E = *nodeIterator;
            int distanceFromEToJ = steinerDistanceTable[E][j];

            // get index of the subset D-E
            cout<<"Getting D-E\n";
            cout << "J: "<< j << endl;
            cout << "Element of D value" << distanceFromEToJ << endl;
            cout << "Value for E is" << E << endl;

            set<int> DMinusE;
            DMinusE = elementOfD;

            for (set<int>::iterator iter = DMinusE.begin(); iter != DMinusE.end();)
            {
              if (*iter == E)
              {
                iter = DMinusE.erase(iter);
              }
              else
              {
                ++iter;
              }
            }
            int dist,indexOfSubsetDMinusE;
            if(DMinusE.size() == 1)
            {
              set<int>::iterator node = DMinusE.begin();
              indexOfSubsetDMinusE = *node;
            }
            else
            {
              indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE);
            }
            // get the real value here
            cout << " indexOfSubsetDMinusE =  " << indexOfSubsetDMinusE <<endl;
            int distanceFromDMinusEToJ = steinerDistanceTable[indexOfSubsetDMinusE][j];
            cout<<"distanceFromDMinusEToJ = "<< distanceFromDMinusEToJ <<"\n";
            dist = distanceFromEToJ + distanceFromDMinusEToJ;
            cout<<"Total Distance = "<< dist <<"\n";
            u = min(u, dist);

            cout << "U is" << u <<endl;
          }

          int val;

          for (int i = 0; i < graph.size(); i++)
          {
            val = graph[i][j];
            if(val != 0 )
            {
              steinerDistanceTable[indexOfSubset][i] = min(steinerDistanceTable[indexOfSubset][i],val + u);
              cout << "I is" << i << endl;
              cout << "Steiner table data is is" << steinerDistanceTable[indexOfSubset][i] << endl;
            }
          }
        }
      }
    }
  }
  int vm = 999;
  int v = 999;

  for (int j = 0; j < graph.size(); j++)
  {
    int u = 999;

    set <int> C1 (C.begin(),C.end());

    for (set<int>::iterator setiter = C1.begin(); setiter != C1.end();setiter++)
    {
      int E = *setiter;
      set<int> subsetCMinusi (C.begin(), C.end());
      int distEJ = graph[E][j];
      for (set<int>::iterator iter = subsetCMinusi.begin(); iter != subsetCMinusi.end();)
      {
        if (*iter == E)
        {
          iter = subsetCMinusi.erase(iter);
        }
        else
        {
          ++iter;
        }
      }
      int indexOfSubsetCMinusi;
      if(subsetCMinusi.size() != 1)
      {
        indexOfSubsetCMinusi = fetchIndexofMapofSets(subsetCMinusi);
      }
      else
      {
        set<int>::iterator node = subsetCMinusi.begin();
        indexOfSubsetCMinusi = *node;
      }

      int val = steinerDistanceTable[indexOfSubsetCMinusi][j];
      cout << " Val: " << val <<endl;
      cout << "E is:  " << E;
      cout << "J: "<<j<<endl;
      cout << " distEJ: " << distEJ <<endl;
      cout << " IndexSubset: " << indexOfSubsetCMinusi;

      u = min(u,distEJ + val );
      
      cout << " U is: " <<u << endl;
    }

  cout << " U final is:" <<u;
  int value = graph[q][j];
  cout<< "value is: " << value <<endl;
  if (value < 999 && value > 0 )
  {
      v = min(v,value+u);
  }
  cout << " V is: " << v << endl;

  }
  for(int i = 0 ; i < steinerDistanceTable.size() ; i++ )
  {
    for(int j = 0; j < graph.size() ; j++)
    {
      cout << " " <<steinerDistanceTable[i][j];
    }
    cout << " " <<endl;
  }


  cout<<"Minimum Steiner Distance"<<v<<"\n";
  return v;
}

void testSubsets()
{
  int n;
  vector<set<int>> allSubsets;
  set<int>::iterator it;

  cout << "Enter size of the set\n";
  cin >> n;

  vector<int> arr(n);

  cout << "Enter Elements of the set\n";
  for (int i = 0; i < n; i++)
    cin >> arr[i];

  // allSubsets = generateSubsets(arr, n);

  for (int i = 0; i < allSubsets.size(); i++)
  {
    set<int> numbersSet = allSubsets[i];
    for (it = numbersSet.begin(); it != numbersSet.end(); it++)
    {
      cout << ' ' << *it;
    }
    cout << "\n";
  }
}

int main()
{
  vector<vector<int> > graph;
  int ySize;

  int a[][10] =  {{0, 3, 1, infi, infi},
                  {3, 0 , 7, 5, 1},
                  {1, 7 , 0, 2, infi},
                  {infi, 5 , 2, 0, 7},
                  {infi,1,infi,7,0}
                  };
  allpairshort(a, 5);

  for (int i = 0; i < 5; i++){        //creating row
    graph.push_back(vector<int>());
  }

  for (int n = 0; n < 5; n++){        //creating columns for the rows
    for (int m = 0; m < 5; m++){
      graph[m].push_back(0);
    }
  }

  for (int m = 0; m < 5; m++){        //storing and printing data
    for (int n = 0; n < 5; n++){
    //  vec[n].push_back(arr[m][n]);
      vec[m][n] = a[m][n];
    }
    cout << "\n";
  }

  cout << "Enter size of terminals\n";
  cin >> ySize;
  vector<int> y(ySize, 0);
  for (int i = 0; i < ySize; i++)
  {
    cin >> y[i];
  }

  int q = y[0];
  y.erase(y.begin());

  set<set<int> > allSubsets = generateSubsets(y, y.size());

  calculateSteiner(y, graph, allSubsets);


  /*for(int i = 0 ; i < 5 ; i++ )
  {
      for(int j = 0; j < 5 ; j++)
      {
          cout << " " <<a[i][j];
      }
      cout << " " <<endl;
  }*/
}
