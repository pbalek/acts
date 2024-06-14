#include "ROOT/RCsvDS.hxx"
#include <string>
#include <unordered_map>

void csv2root(std::string fileCSV="itk_output/event000000000-seed.csv", std::string fileROOT="event000000000-seed.root")
{
    // auto rdf = ROOT::RDF::FromCSV(fileCSV, true, ',', -1, {{"measurement_id_1",'L'},{"measurement_id_2",'L'},
    //     {"measurement_id_3",'L'},{"geometry_id_1",'L'},{"geometry_id_2",'L'},{"geometry_id_3",'L'},
    //     {"x_1",'D'},{"x_2",'D'},{"x_3",'D'},
    //     {"y_1",'D'},{"y_2",'D'},{"y_3",'D'},
    //     {"z_1",'D'},{"z_2",'D'},{"z_3",'D'},
    //     {"var_r_1",'D'},{"var_r_2",'D'},{"var_r_3",'D'},
    //     {"var_z_1",'D'},{"var_z_2",'D'},{"var_z_3",'D'},
    //     {"z_vertex",'D'},{"seed_quality",'D'}});

    // auto rdf = ROOT::RDF::FromCSV(fileCSV, true, ',', -1, {{"z_3",'D'},{"z_2",'D'}});

    auto rdf = ROOT::RDF::FromCSV(fileCSV, true, ',', -1);

    if(rdf.HasColumn("z_3"))
    {
        std::cout<<"column z_3 present"<<std::endl;
    }
    else
    {
        std::cout<<"column z_3 missing"<<std::endl;
    }
    
    rdf.Snapshot("tree", fileROOT);
    return;
}


