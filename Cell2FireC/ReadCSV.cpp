#include "ReadCSV.h"
#include "fuelmodelBP.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#include <map> // MHK
#include <sstream> // For K opt

/*
 * Constructur
 */
CSVReader::CSVReader(std::string filename, std::string delm){
	this->fileName = filename;
	this->delimeter = delm;		
}
 
 
/*
* Parses through csv file line by line and returns the data
* in vector of vector of strings.
*/
std::vector<std::vector<std::string>> CSVReader::getData(){
	std::ifstream file(this->fileName);
	std::vector<std::vector<std::string> > dataList;
	std::string line = "";
	
	// Iterate through each line and split the content using delimeter
	while (getline(file, line))
	{
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(this->delimeter));
		dataList.push_back(vec);
	}
	// Close the File
	file.close();
 
	return dataList;
}


/*
* Prints data to screen inside the DF obtained from the CSV file
*/
void CSVReader::printData(std::vector<std::vector<std::string>> & DF){
	// Print the content of row by row on screen
	for(std::vector<std::string> vec : DF)
	{
		for(std::string data : vec)
		{
			if(data.length() >0) std::cout << " " << data << " ";
			else std::cout << " " <<  "NaN" << " ";
		}
		std::cout<<std::endl;
	}
}


/*
* Populates the df input objects based on the DF csv file for each row/cell (spanish version)
*/
void CSVReader::parseDF(inputs * df_ptr, std::vector<std::vector<std::string>> & DF, int NCells){
	int i;
	
	// Floats 
	float cur, elev, ws, waz, saz, cbd, cbh, ccf, ps;
	
	// Integers 
	int nftype;
	
	// CChar
	const char * faux;
	std::string::size_type sz;   // alias of size_t
	
	// Loop over cells (populating per row)
	for (i=1; i <= NCells; i++){
		//printf("Populating DF for cell %d\n", i);
		faux = DF[i][0].append(" ").c_str();
		
		if (DF[i][3].compare("") == 0) elev = 0;
		else elev = std::stof (DF[i][3] ,&sz);
		
		if (DF[i][4].compare("") == 0) ws = 0;
		else ws = std::stof (DF[i][4], &sz);
		
		if (DF[i][5].compare("") == 0) waz = 0;
		else waz = std::stoi (DF[i][5] ,&sz) + 180.;// + 2*90;  // CHECK!!!!
		if (waz >= 360) waz = waz - 360;
		
		if (DF[i][6].compare("") == 0) ps = 0;
		else ps = std::stof (DF[i][6] ,&sz);
		
		if (DF[i][7].compare("") == 0) saz = 0;
		else saz = std::stof (DF[i][7] ,&sz);
		
		if (DF[i][8].compare("") == 0) cur = 0;
		else cur = std::stof (DF[i][8] ,&sz);
			
		if (DF[i][9].compare("") == 0) cbd = 0;
		else cbd = std::stof (DF[i][9], &sz);
		
		if (DF[i][10].compare("") == 0) cbh = 0;
		else cbh = std::stof (DF[i][10], &sz);
		
		if (DF[i][11].compare("") == 0) ccf = 0;
		else ccf = std::stof (DF[i][11], &sz);
		
		if (DF[i][12].compare("") == 0) nftype = 0;
		else nftype = std::stoi (DF[i][12], &sz);
		
		// Set values
		strncpy(df_ptr->fueltype, faux, 4);
		df_ptr->elev=elev; df_ptr->ws=ws; df_ptr->waz=waz;  
		df_ptr->ps=ps; df_ptr->saz=saz; df_ptr->cur=cur; 
		df_ptr->cbd=cbd;df_ptr->cbh=cbh; df_ptr->ccf=ccf;
		df_ptr->nftype=nftype;
			
		// Next pointer
		df_ptr++;
	}
}


/*
* Populates vector of size NCells with type number based on lookup table (Spain version)
*/
void CSVReader::parseNDF(std::vector<int> &NFTypes, std::vector<std::vector<std::string>> & DF, int NCells){
	int i;
	
	// Ints 
	int FType;
	
	// CChar
	const char * faux;
	std::string::size_type sz;   // alias of size_t
	
	// Loop over cells (populating per row)
	for (i=1; i <= NCells; i++){
		//printf("Populating DF for cell %d\n", i);
		if (DF[i][12].compare("") == 0) FType = 0;
		else FType = std::stoi (DF[i][12], &sz);
			
		// Set values
		NFTypes.push_back(FType);
		
	}
}


/*
* Populate Weather DF Spain
*/
void CSVReader::parseWeatherDF(weatherDF * wdf_ptr, std::vector<std::vector<std::string>> & DF, int WPeriods){
	int i;
	
	//Strings
	std::string instance, datetime;
	std::string::size_type sz;   // alias of size_t
	
	//Floats 
	float ws, waz, slope;
	
	//Floats for MC 
	float mc1, mc10, mc100, mcWoody, mcHerb;
	
	// Loop over cells (populating per row)
	for (i=1; i <= WPeriods; i++){
		//printf("Populating Weather DF period %d\n", i);
		instance = DF[i][0];
		datetime = DF[i][1];
		

		if (DF[i][3].compare("") == 0) waz = 0;
		else {waz = std::stoi (DF[i][3] ,&sz); //+ 180/2;   // DEBUGGING THE ANGLE 
			if (waz >= 360){
				waz = waz - 360;
			}
		}
		
		// WS [2]
		if (DF[i][2].compare("") == 0) ws = 0;
		else ws = std::stof (DF[i][2], &sz);
		
		slope = std::stof (DF[i][4], &sz);
		mc1 = std::stof (DF[i][5], &sz);
		mc10 = std::stof (DF[i][6], &sz);
		mc100 = std::stof (DF[i][7], &sz);
		mcHerb = std::stof (DF[i][8], &sz);
		mcWoody = std::stof (DF[i][9], &sz);


		// // Slope [4]
		// if (DF[i][4].compare("") == 0) slope = 0;
		// else slope = std::stof (DF[i][4], &sz);
		
		// // MC1 [5]
		// if (DF[i][5].compare("") == 0) mc1 = 0;
		// else mc1 = std::stoi (DF[i][5], &sz);

		// // MC10 [6]
		// if (DF[i][6].compare("") == 0) mc10 = 0;
		// else mc10 = std::stoi (DF[i][6], &sz);

		// // MC100 [7]
		// if (DF[i][7].compare("") == 0) mc100 = 0;
		// else mc100 = std::stoi (DF[i][7], &sz);

		// // MC Woody [8]
		// if (DF[i][8].compare("") == 0) mcWoody = 0;
		// else mcWoody = std::stoi (DF[i][8], &sz);

		// // MC Herb [9]
		// if (DF[i][8].compare("") == 0) mcHerb = 0;
		// else mcHerb = std::stoi (DF[i][9], &sz);

		// if (DF[i][4].compare("") == 0) scenario = 0;
		// else scenario = std::stof (DF[i][4], &sz);
		
		// Set values
		wdf_ptr->instance = instance;
		wdf_ptr->datetime = datetime;
		wdf_ptr->ws=ws; 
		wdf_ptr->waz=waz; 
		wdf_ptr->slope=slope;
		wdf_ptr->mc1=mc1;
		wdf_ptr->mc10=mc10;
		wdf_ptr->mc100=mc100;
		wdf_ptr->mcWoody=mcWoody;
		wdf_ptr->mcHerb=mcHerb;
		// wdf_ptr->scenario = scenario;
			
		// Next pointer
		wdf_ptr++;
	}
	
}


/*
* Populate K opt (MHK)
*/
// void CSVReader::parseKopt(koptDF * kopt_ptr, std::vector<std::vector<std::string>> & DF){
//     std::string filename = "k_opt.csv"; // Set filename to be from args
//     std::ifstream file(filename);

//     std::vector<float> var1, var2, var3, var4, var5, var6, var7;

//     // Read and discard the first line (Column names)
//     std::string line;
//     getline(file, line);

//     // if (file.is_open()) {
//     while (getline(file, line)) {
//         std::stringstream ss(line);
//         std::string value;

//         // Read in each value and store it in the appropriate vector
//         getline(ss, value, ',');
//         var1.push_back(stof(value));
//         getline(ss, value, ',');
//         var2.push_back(stof(value));
//         getline(ss, value, ',');
//         var3.push_back(stof(value));
//         getline(ss, value, ',');
//         var4.push_back(stof(value));
//         getline(ss, value, ',');
//         var5.push_back(stof(value));
//         getline(ss, value, ',');
//         var6.push_back(stof(value));
//         getline(ss, value, ',');
//         var7.push_back(stof(value));
    
// 		// Set values
// 		kopt_ptr->ecc=var1;
// 		kopt_ptr->K1=var3;
// 		kopt_ptr->K2=var4; 
// 		kopt_ptr->K3=var5; 
// 		kopt_ptr->K4=var6;
// 		kopt_ptr->K5=var7;
			
// 		// Next pointer
// 		kopt_ptr++;
//     }

// }

/*
* Populate BBO Tuning factors
*/
// void CSVReader::parseKoptDict2(std::unordered_map<int, std::vector<float>> & bbo, std::vector<std::vector<std::string>> & DF, int NFTypes){
// 	// Integers
// 	int i, j, ftype;
// 	int ffactors = 4;
// 	std::vector<float> bboFactors;
// 	std::string::size_type sz;   // alias of size_t
	
// 	// Loop over cells (populating per row)
// 	for (i=1; i <= NFTypes; i++){
// 		// Clean the vector before the fuels
// 		bboFactors.clear();
		
// 		//DEBUGprintf("Populating Ignition points: %d\n", i);
// 		ftype = std::stoi(DF[i][0], &sz);
		
// 		for (j=1; j <= ffactors; j++){
// 			bboFactors.push_back(std::stof(DF[i][j], &sz));
// 		}

// 		//Set values							 
// 		bbo.insert(std::make_pair(ftype, bboFactors));	
// 	}
	
	
// }

void CSVReader::parseKoptDict(std::string filename, std::unordered_map<double, std::vector<double>> *kdict_ptr) {
    if (!kdict_ptr) {
        std::cerr << "Error: The input pointer is null." << std::endl;
        return;
    }
    
    // Open the CSV file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return;
    }

    // Read and discard the first line (Column names)
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        // Split the line into values
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }

        // Add the row to the dictionary if it's not empty
        // std::cout << "DEBUGGING ReadCSV.cpp Adding rows to dict" << std::endl;
        if (!row.empty()) {
            double key = 0.0; // Initialize key to a default value
            if (!row.empty()) {
                key = row[0];
                // std::cout << "DEBUGGING ReadCSV.cpp P1" << std::endl;
                row.erase(row.begin()); // Remove the key from the row
                // std::cout << "DEBUGGING ReadCSV.cpp P2" << std::endl;
            }
            // std::cout << "DEBUGGING ReadCSV.cpp Finish P2" << std::endl;

            (*kdict_ptr)[key] = row;

            // std::cout << "DEBUGGING ReadCSV.cpp P3" << std::endl;
        }
    }

    file.close();

	// std::cout << "DEBUGGING ReadCSV.cpp Read out file" << std::endl;
    
	// Print the dictionary values
    // std::cout << "Parsed Dictionary:\n";
    // for (auto iter = kdict_ptr->begin(); iter != kdict_ptr->end(); ++iter) {
    //     std::cout << iter->first << ": [";
    //     for (double d : iter->second) {
    //         std::cout << d << ", ";
    //     }
    //     std::cout << "]\n";
    // }
}

    // // Free the dictionary memory
    // delete myDictionary;


// 	std::ifstream file(filename);
//     std::string line;

// 	std::cout << "DEBUGGING ReadCSV.cpp Read in file" << std::endl;

//     // Read and discard the first line (Column names)
//     getline(file, line);
	
// 	std::cout << "DEBUGGING ReadCSV.cpp Reading" << std::endl;

//     if (file.is_open()) {
//         while (std::getline(file, line)) {
//             std::istringstream iss(line);
//             double key;
//             std::vector<double> values;
//             if (!(iss >> key)) {
//                 continue;
//             }
//             double value;
//             while (iss >> value) {
//                 values.push_back(value);
//                 if (iss.peek() == ',') {
//                     iss.ignore();
//                 }
//             }
//             (*koptDict)[key] = values;
//         }
//         file.close();
// 	std::cout << "DEBUGGING ReadCSV.cpp Done reading" << std::endl;
//     }
// }


/*
* Populate IgnitionDF
*/
void CSVReader::parseIgnitionDF(std::vector<int> & ig, std::vector<std::vector<std::string>> & DF, int IgPeriods){
	// Integers
	int i, igcell;
	std::string::size_type sz;   // alias of size_t
	
	// Loop over cells (populating per row)
	for (i=1; i <= IgPeriods; i++){
		//DEBUGprintf("Populating Ignition points: %d\n", i);
		igcell = std::stoi(DF[i][1], &sz);
		
		// Set values
		ig[i-1]= igcell;
			
		// Next pointer
		//ig_ptr++;
	}
	
	
}


/*
* Populate HarvestedDF
*/
void CSVReader::parseHarvestedDF(std::unordered_map<int, std::vector<int>> & hc, std::vector<std::vector<std::string>> & DF, int HPeriods){
	// Integers
	int i, j, hcell;
	std::vector<int> toHarvestCells;
	std::string::size_type sz;   // alias of size_t
	
	// Loop over cells (populating per row)
	for (i=1; i <= HPeriods; i++){
		// Clean the vector before the new year 
		toHarvestCells.clear();
		
		// Loop over years of the simulation
		for(j=1; j < DF[i].size(); j++){
			hcell = std::stoi(DF[i][j], &sz);
		
			// Set values
			toHarvestCells.push_back(hcell);
		}
		
		// Populate unordered set 
		hc.insert(std::make_pair(i, toHarvestCells));	
		
	}
	
	
}


/*
* Populate BBO Tuning factors
*/
void CSVReader::parseBBODF(std::unordered_map<int, std::vector<float>> & bbo, std::vector<std::vector<std::string>> & DF, int NFTypes){
	// Integers
	int i, j, ftype;
	int ffactors = 4;
	std::vector<float> bboFactors;
	std::string::size_type sz;   // alias of size_t
	
	// Loop over cells (populating per row)
	for (i=1; i <= NFTypes; i++){
		// Clean the vector before the fuels
		bboFactors.clear();
		
		//DEBUGprintf("Populating Ignition points: %d\n", i);
		ftype = std::stoi(DF[i][0], &sz);
		
		for (j=1; j <= ffactors; j++){
			bboFactors.push_back(std::stof(DF[i][j], &sz));
		}

		//Set values							 
		bbo.insert(std::make_pair(ftype, bboFactors));	
	}
	
	
}


void CSVReader::parseForestDF(forestDF * frt_ptr, std::vector<std::vector<std::string>> & DF){
	// Ints 
	int cellside, rows, cols;
	int i, j;
	
	// Others 
	std::vector<std::unordered_map<std::string, int>> adjCells;
 	std::string::size_type sz;   // alias of size_t
	std::vector<std::vector<int>> coordCells;
	std::unordered_map<std::string, int> Aux;
	std::vector<int> Aux2;
	
	std::string North = "N";
    std::string South = "S";
    std::string East = "E";
    std::string West = "W";
    std::string NorthEast = "NE";
    std::string NorthWest = "NW";
    std::string SouthEast = "SE";
    std::string SouthWest = "SW";
	
	// Filling DF
	//DEBUGprintf("Populating Forest DF\n");
	
	cols = std::stoi (DF[0][1], &sz);
	rows = std::stoi (DF[1][1], &sz);
	cellside = std::stoi (DF[4][1], &sz);
	
	//DEBUGprintf("cols: %d,  rows:  %d,   cellside:  %d\n", cols, rows, cellside);
	
	// CoordCells and Adjacents
	int n = 1; 
	int r, c;
	//std::cout  << "CoordCells Debug" << std::endl;
	for (r=0; r<rows; r++){
		for (c=0; c < cols; c++){
			
			/*   CoordCells  */
			Aux2 = std::vector<int>();
			Aux2.push_back(c); 
            Aux2.push_back(rows-r-1);   
			coordCells.push_back(Aux2);                    
			//printf("i,j = %d,%d\n", r,c);
			//std::cout << "x: " << coordCells[c + r*(cols)][0] <<  "  y: " << coordCells[c + r*(cols)][1]  <<   std::endl;
					
			/*   Adjacents  */
			// if we have rows (not a forest = line)
			if (rows>1){
				
				// Initial row
				if(r == 0){
					
					if (c == 0){
                        Aux = {{North,-1},{NorthEast,-1},{NorthWest,-1},{South,n+cols},{SouthEast,n+cols+1}, 
							        {SouthWest,-1}, {East,n+1},{West,-1}};
                        adjCells.push_back(Aux);
						n++;
					}
                    if (c == cols - 1){
                        Aux = {{North,-1},{NorthEast,-1},{NorthWest,-1},{South, n+cols},{SouthEast,-1},
										{SouthWest, n+cols-1,}, {East,-1}, {West,n-1}};
						adjCells.push_back(Aux);
                        n++;
					}
                    if (c > 0 && c < cols-1){    
                        Aux = {{North, -1},{NorthEast,-1},{NorthWest,-1},{South,n+cols},{SouthEast,n+cols+1}, 
									{SouthWest, n+cols-1}, {East, n+1},{West,n-1}};
						adjCells.push_back(Aux);
						n++;
					}
				}
				
				// In between
				if (r > 0 && r < rows - 1){
                    if (c == 0){
                        Aux = {{North, n-cols} , {NorthEast, n-cols+1 }, {NorthWest,-1}, {South, n+cols}, 
									{SouthEast, n+cols+1} , {SouthWest,-1}, {East, n+1} ,{West,-1}};
						adjCells.push_back(Aux);
                        n++;
					}
                    if (c == cols-1){
                        Aux = {{North, n-cols}, {NorthEast,-1}, {NorthWest, n-cols-1},{South, n+cols}, 
									{SouthEast,-1}, {SouthWest, n+cols-1}, {East,-1}, {West, n-1}};
                        adjCells.push_back(Aux);
						n++;
					}
                    if (c>0 && c<cols-1){    
                        Aux = {{North, n-cols}, {NorthEast, n-cols+1} , {NorthWest, n-cols-1}, {South, n+cols}, 
									{SouthEast, n+cols+1} , {SouthWest, n+cols-1}, {East, n+1}, {West, n-1}};
						adjCells.push_back(Aux);
                        n++;    
					}
				}
				
				// Final row
				if (r == rows-1){
                    if (c == 0){
                        Aux = {{North,n-cols}, {NorthEast,n-cols+1}, {NorthWest,-1}, {South,-1}, {SouthEast,-1}, 
									{SouthWest,-1,}, {East,n+1}, {West,-1}};
						adjCells.push_back(Aux);				 
                        n++;    
					}
                        
                    if (c == cols-1){
                        Aux = {{North,n-cols}, {NorthEast,-1}, {NorthWest,n-cols-1}, {South,-1}, {SouthEast,-1}, 
									{SouthWest,-1}, {East,-1}, {West,n-1}};
						adjCells.push_back(Aux);
                        n++;    
					}
                    if (c>0 and c<cols-1){    
                        Aux = {{North,n-cols}, {NorthEast, n-cols+1}, {NorthWest,n-cols-1}, {South,-1}, 
									{SouthEast,-1} , {SouthWest,-1}, {East,n+1}, {West,n-1}};
						adjCells.push_back(Aux);
						n++;    
					}
				
				}
			}	
				
			// One line
			if (rows == 1){
				if (c == 0){
					Aux = {{North,-1}, {NorthEast,-1}, {NorthWest,-1}, {South,-1}, {SouthEast,-1}, 
								{SouthWest,-1}, {East,n+1}, {West,-1}};
					adjCells.push_back(Aux);
					n++;    
				}
				if (c == cols-1){
					Aux = {{North,-1}, {NorthEast,-1}, {NorthWest,-1}, {South,-1}, {SouthEast,-1}, 
								{SouthWest,-1}, {East,-1},{West,n-1}};
					adjCells.push_back(Aux);
					n++;    
				}
				if (c>0 && c<cols-1){						
					Aux = {{North,-1}, {NorthEast,-1}, {NorthWest,-1}, {South,-1}, {SouthEast,-1}, 
								{SouthWest,-1}, {East,n+1}, {West,n-1}};
					adjCells.push_back(Aux);
					n++;    
				}
			}
		}
	}
	
	
	
	// Adjacents cells
	//std::cout  << "Adjacents Debug" << std::endl;
	/*for (i=0; i<adjCells.size();i++){
		std::cout << "Cell "<< i+1 << " =  "; 
		for (auto & nb : adjCells[i]){
			std::cout << " " << nb.first << " : " << nb.second;
		}
		std::cout << std::endl;
	}
	*/
	
	
	// Set values
	frt_ptr->cellside = cellside;
	frt_ptr->rows = rows;
	frt_ptr->cols = cols;
	frt_ptr->coordCells = coordCells;
	frt_ptr->adjCells = adjCells;
		
		
	
}


void CSVReader::printDF(inputs df){
	std::cout << df.fueltype; std::cout << " ";
	std::cout << " " << df.elev; std::cout << " " << df.ws; std::cout << " " << df.waz; 
	std::cout << " " << df.ps; std::cout << " " << df.saz; std::cout << " " << df.cur; 
	std::cout << " " << df.cbd; std::cout << " " << df.cbh; std::cout << " " << df.ccf << std::endl;
}


void CSVReader::printWeatherDF(weatherDF wdf){
	std::cout << " " << wdf.datetime; 
	std::cout << " " << wdf.ws; std::cout << " " << wdf.waz; 
}

/*
int main()
{
	// Creating an object of CSVWriter
	CSVReader reader("example.csv");
 
	// Get the data from CSV File
	std::vector<std::vector<std::string> > dataList = reader.getData();
 
	// Print the content of row by row on screen
	for(std::vector<std::string> vec : dataList)
	{
		for(std::string data : vec)
		{
			std::cout<<data<< " ";
		}
		std::cout<<std::endl;
	}
	return 0;
 
}
*/
