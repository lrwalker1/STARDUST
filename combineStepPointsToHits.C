#include <cmath>
#include <random>

// Function to apply Gaussian blur
double applyGaussianBlur(double energy, double fwhm) {
    static std::default_random_engine generator;
    double sigma = fwhm / 2.35482; // Convert FWHM to sigma
    std::normal_distribution<double> distribution(energy, sigma);
    return distribution(generator);
}

// ROOT script to combine nearby Geant4 step points to hits in a detector
// based on the space resolution of the detector
void combineStepPointsToHits(
        const char *input="output2.root",
        const char *output="hits2.root",
        double resolution=0.78125/*mm*/)
{
    // input, which is GEARS output
    TChain *t = new TChain("t");
    t->Add(input);
    int nStepPoints; // number of Geant4 step points
    t->SetBranchAddress("n",&nStepPoints);
    // parameters of step points
    vector<double> *xx=0, *yy=0, *zz=0, *de=0, *time=0, *k=0;
    vector<int> *vlm=0, *pdg=0, *pro=0; // copy number of a Geant4 volume, particle ID, and process ID
    TBranch *bx, *by, *bz, *be, *bv, *bt, *bpdg, *bpro, *bk;
    t->SetBranchAddress("xx",&xx, &bx); // local x
    t->SetBranchAddress("yy",&yy, &by); // local y
    t->SetBranchAddress("zz",&zz, &bz); // local z
    t->SetBranchAddress("de",&de, &be); // energy deposition
    t->SetBranchAddress("vlm",&vlm, &bv);
    t->SetBranchAddress("t",&time, &bt); // time
    t->SetBranchAddress("pdg",&pdg,&bpdg); // particle ID
    t->SetBranchAddress("pro",&pro,&bpro); // process ID
    t->SetBranchAddress("k",&k,&bk); // kinetic energy

    // output, which contains a tree filled with combined hits
    TFile *file = new TFile(output, "recreate");
    TTree *tree = new TTree("h","combined hits");
    int n; // number of combined hits
    int evt; // id of event from Geant4 simulation
    tree->Branch("n",  &n,  "n/I");
    tree->Branch("evt",&evt,"evt/I");
    // parameters of combined hits
    double x[1000], y[1000], z[1000], e[1000], times[1000], edet[1000];
    int pdg_combined[1000], vlm_combined[1000];
    tree->Branch("x",x,"x[n]/D");
    tree->Branch("y",y,"y[n]/D");
    tree->Branch("z",z,"z[n]/D");
    tree->Branch("e",e,"e[n]/D");
    tree->Branch("t",times,"t[n]/D");
    tree->Branch("pdg",pdg_combined,"pdg[n]/I");
    tree->Branch("vlm",vlm_combined,"vlm[n]/I");
    tree->Branch("edet",edet,"edet[n]/D");

    // File to write verbose information
    // ofstream debugFile("19.txt");

    // main loop to combine step points
    double dx=0, dy=0, dz=0, dr=0; // distances between step points
    int nevt = t->GetEntries(); // total number of events simulated
    cout<<nevt<<" events to be processed"<<endl;
    for (evt=0; evt<nevt; evt++) {
        if (evt%10000==0) cout<<evt<<" events processed"<<endl;
        t->GetEntry(evt); // get information of step points from input tree

        n = 0; // reset counter of combined hits for a new event
        // vector<int> debugSteps; // reset debug steps for a new event
        vector<int> contributingSteps; // steps contributing to the current hit

        for (int i=0; i<nStepPoints; i++) { // loop over step points
            if (de->at(i)==0) continue; // skip step points with no energy deposition
            if (vlm->at(i)<1) continue; // skip step points with vlm less than 1
            
            // Check if the current step point belongs to HPGe or DSSD detectors
            bool isHPGe = (vlm->at(i) >= 1 && vlm->at(i) <= 21);
            bool isDSSD = (vlm->at(i) == 801 || vlm->at(i) == 802);

            // debugSteps.push_back(i); // add step index to debugSteps

            if (n == 0) { // no combined hit yet, create the 1st one
                x[n] = xx->at(i); y[n] = yy->at(i); z[n] = zz->at(i); e[n] = de->at(i); times[n] = time->at(i); pdg_combined[n] = pdg->at(i); vlm_combined[n] = vlm->at(i);
                contributingSteps.push_back(i); // add step index to contributingSteps
                n++; // increase the hit index by 1
                continue;
            }

            if (isHPGe) {
                // HPGe Detectors (1-21): Sum energy, use initial time, first PDG, and average position
                if (vlm_combined[n-1] == vlm->at(i)) { // same HPGe detector
                    // get energy weighted position
                    x[n-1] = (x[n-1] * e[n-1] + xx->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    y[n-1] = (y[n-1] * e[n-1] + yy->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    z[n-1] = (z[n-1] * e[n-1] + zz->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    e[n-1] += de->at(i); // deposited energies need to be summed
                    contributingSteps.push_back(i); // add step index to contributingSteps
                } else { // create a new hit
                    // Check if the previous hit meets the energy and detector criteria before creating a new hit
                    // if (e[n-1] >= 18 && e[n-1] <= 21 && vlm_combined[n-1] >= 1 && vlm_combined[n-1] <= 21) {
                    //     debugFile << "Event: " << evt << "\n";
                    //     debugFile << "Step points contributing to hit:\n";
                    //     for (auto stepIdx : contributingSteps) {
                    //         debugFile << "Step " << stepIdx << ": "
                    //                   << "xx=" << xx->at(stepIdx) << ", "
                    //                   << "yy=" << yy->at(stepIdx) << ", "
                    //                   << "zz=" << zz->at(stepIdx) << ", "
                    //                   << "de=" << de->at(stepIdx) << ", "
                    //                   << "t=" << time->at(stepIdx) << ", "
                    //                   << "vlm=" << vlm->at(stepIdx) << ", "
                    //                   << "pdg=" << pdg->at(stepIdx) << ", "
                    //                   << "pro=" << pro->at(stepIdx) << ", "
                    //                   << "k=" << k->at(stepIdx) << "\n";
                    //     }
                    //     debugFile << "\nPrevious steps:\n";
                    //     for (int j = 0; j < contributingSteps.front(); j++) {
                    //         debugFile << "Step " << j << ": "
                    //                   << "xx=" << xx->at(j) << ", "
                    //                   << "yy=" << yy->at(j) << ", "
                    //                   << "zz=" << zz->at(j) << ", "
                    //                   << "de=" << de->at(j) << ", "
                    //                   << "t=" << time->at(j) << ", "
                    //                   << "vlm=" << vlm->at(j) << ", "
                    //                   << "pdg=" << pdg->at(j) << ", "
                    //                   << "pro=" << pro->at(j) << ", "
                    //                   << "k=" << k->at(j) << "\n";
                    //     }
                    //     debugFile << "\n";
                    // }
                    // Reset contributingSteps for the new hit
                    contributingSteps.clear();
                    x[n] = xx->at(i); y[n] = yy->at(i); z[n] = zz->at(i); e[n] = de->at(i); times[n] = time->at(i); pdg_combined[n] = pdg->at(i); vlm_combined[n] = vlm->at(i);
                    n++;
                    contributingSteps.push_back(i); // add step index to contributingSteps
                }
            } else if (isDSSD) {
                // DSSD Detectors (801-802): Sum energy within resolution distance, use initial time, first PDG, and average position
                dx = xx->at(i) - x[n-1]; dy = yy->at(i) - y[n-1]; dz = zz->at(i) - z[n-1];
                dr = sqrt(dx*dx + dy*dy + dz*dz);
                if (dr > resolution) { // create a new hit far away from the previous one
                    // Check if the previous hit meets the energy and detector criteria before creating a new hit
                    // if (e[n-1] >= 18 && e[n-1] <= 21 && vlm_combined[n-1] >= 1 && vlm_combined[n-1] <= 21) {
                    //     debugFile << "Event: " << evt << "\n";
                    //     debugFile << "Step points contributing to hit:\n";
                    //     for (auto stepIdx : contributingSteps) {
                    //         debugFile << "Step " << stepIdx << ": "
                    //                   << "xx=" << xx->at(stepIdx) << ", "
                    //                   << "yy=" << yy->at(stepIdx) << ", "
                    //                   << "zz=" << zz->at(stepIdx) << ", "
                    //                   << "de=" << de->at(stepIdx) << ", "
                    //                   << "t=" << time->at(stepIdx) << ", "
                    //                   << "vlm=" << vlm->at(stepIdx) << ", "
                    //                   << "pdg=" << pdg->at(stepIdx) << ", "
                    //                   << "pro=" << pro->at(stepIdx) << ", "
                    //                   << "k=" << k->at(stepIdx) << "\n";
                    //     }
                    //     debugFile << "\nPrevious steps:\n";
                    //     for (int j = 0; j < contributingSteps.front(); j++) {
                    //         debugFile << "Step " << j << ": "
                    //                   << "xx=" << xx->at(j) << ", "
                    //                   << "yy=" << yy->at(j) << ", "
                    //                   << "zz=" << zz->at(j) << ", "
                    //                   << "de=" << de->at(j) << ", "
                    //                   << "t=" << time->at(j) << ", "
                    //                   << "vlm=" << vlm->at(j) << ", "
                    //                   << "pdg=" << pdg->at(j) << ", "
                    //                   << "pro=" << pro->at(j) << ", "
                    //                   << "k=" << k->at(j) << "\n";
                    //     }
                    //     debugFile << "\n";
                    // }
                    // Reset contributingSteps for the new hit
                    contributingSteps.clear();
                    x[n] = xx->at(i); y[n] = yy->at(i); z[n] = zz->at(i); e[n] = de->at(i); times[n] = time->at(i); pdg_combined[n] = pdg->at(i); vlm_combined[n] = vlm->at(i);
                    n++;
                    contributingSteps.push_back(i); // add step index to contributingSteps
                } else { // combine a nearby step point with the previously combined hit
                    // get energy weighted position
                    x[n-1] = (x[n-1] * e[n-1] + xx->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    y[n-1] = (y[n-1] * e[n-1] + yy->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    z[n-1] = (z[n-1] * e[n-1] + zz->at(i) * de->at(i)) / (e[n-1] + de->at(i));
                    e[n-1] += de->at(i); // deposited energies need to be summed
                    contributingSteps.push_back(i); // add step index to contributingSteps
                }
            }

            // Keep the first time and first PDG
            times[n-1] = times[n-1]; // already set in the first step
            pdg_combined[n-1] = pdg_combined[n-1]; // already set in the first step
            vlm_combined[n-1] = vlm_combined[n-1]; // already set in the first step
        }

        // Check the last hit if it meets the energy and detector criteria
        // if (e[n-1] >= 18 && e[n-1] <= 21 && vlm_combined[n-1] >= 1 && vlm_combined[n-1] <= 21) {
        //     debugFile << "Event: " << evt << "\n";
        //     debugFile << "Step points contributing to hit:\n";
        //     for (auto stepIdx : contributingSteps) {
        //         debugFile << "Step " << stepIdx << ": "
        //                   << "xx=" << xx->at(stepIdx) << ", "
        //                   << "yy=" << yy->at(stepIdx) << ", "
        //                   << "zz=" << zz->at(stepIdx) << ", "
        //                   << "de=" << de->at(stepIdx) << ", "
        //                   << "t=" << time->at(stepIdx) << ", "
        //                   << "vlm=" << vlm->at(stepIdx) << ", "
        //                   << "pdg=" << pdg->at(stepIdx) << ", "
        //                   << "pro=" << pro->at(stepIdx) << ", "
        //                   << "k=" << k->at(stepIdx) << "\n";
        //     }
        //     debugFile << "\nPrevious steps:\n";
        //     for (int j = 0; j < contributingSteps.front(); j++) {
        //         debugFile << "Step " << j << ": "
        //                   << "xx=" << xx->at(j) << ", "
        //                   << "yy=" << yy->at(j) << ", "
        //                   << "zz=" << zz->at(j) << ", "
        //                   << "de=" << de->at(j) << ", "
        //                   << "t=" << time->at(j) << ", "
        //                   << "vlm=" << vlm->at(j) << ", "
        //                   << "pdg=" << pdg->at(j) << ", "
        //                   << "pro=" << pro->at(j) << ", "
        //                   << "k=" << k->at(j) << "\n";
        //     }
        //     debugFile << "\n";
        // }

        if (n > 0) {
            // Apply Gaussian blur to the energy values based on detector type
            for (int i = 0; i < n; i++) {
                if (vlm_combined[i] >= 1 && vlm_combined[i] <= 21) {
                    edet[i] = applyGaussianBlur(e[i], 0.2 * e[i] / 100.0);
                } else if (vlm_combined[i] == 801 || vlm_combined[i] == 802) {
                    edet[i] = applyGaussianBlur(e[i], 0.4 * e[i] / 100.0);
                } else {
                    edet[i] = e[i]; // No blur for other detectors
                }
            }
            tree->Fill(); // insert x, y, z, e, t, pdg values to the output tree
        }
    }

    // save the output tree
    tree->Write("", TObject::kWriteDelete); // write tree, then delete previous
    file->Close(); // close output file
    // debugFile.close(); // close debug file
}

