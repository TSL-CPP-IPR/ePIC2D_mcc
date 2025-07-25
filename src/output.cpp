#include "output.h"


Output::Output(const std::filesystem::path& outputfolder, Domain& domain) : outputfolder(outputfolder), domain(domain) 
{

    std::filesystem::remove_all(outputfolder);  // Clear previous output
    std::filesystem::create_directories(outputfolder);
    //std::filesystem::create_directories(outputfolder / "files");

    //file_data.open(outputfolder / "files" / "results_" + std::to_string(domain.alpha) + ".txt");
    file = H5File(outputfolder / "result.h5", H5F_ACC_TRUNC);
    
    if (file.getId() < 0) 
    {
        throw std::runtime_error("Error opening HDF5 file");
    }
    //create groups
    //particle_group1 = file.createGroup("/electron");
    //particle_group2 = file.createGroup("/ion");
    field_data_group = file.createGroup("/fielddata");
    time_group = file.createGroup("/time_var");
    metadata_group = file.createGroup("/metadata");
    metadata_species = file.createGroup("/metadata_species");

    int sp_no = domain.species_no;
    int t_step = int(domain.NUM_TS/domain.write_interval) + 1;
    
    store_ke = Matrix<double>(t_step,3*sp_no + 3); //1 time coloumn + two field energy coulumn
    store_m = Matrix<double>(t_step,sp_no + 1);
}

void Output::write_particle_data(int ts, Species& species)
{
    std::string group_name = "particle_" + species.name;
    Group particle_group;

    // Check if the group already exists in the map
    auto it = particle_groups.find(species.name);
    if (it != particle_groups.end()) 
    {
        // Group already exists, retrieve it from the map
        particle_group = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        particle_group = file.createGroup(group_name);
        particle_groups[species.name] = particle_group;
    }

    // Create datasets for particle positions and velocities
    hsize_t dims_pos[2] = {species.part_list.size(), 2};
    hsize_t dims_vel[2] = {species.part_list.size(), 2};
    hsize_t Rank = 2;

    DataSpace dataspace_pos(Rank, dims_pos);
    DataSpace dataspace_vel(Rank, dims_vel);

    H5::DataSet dataset_pos = particle_group.createDataSet("pos" + std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_pos);
    H5::DataSet dataset_vel = particle_group.createDataSet("vel" + std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_vel);

    // Write particle positions and velocities to the datasets
    std::vector<double> positions;
    std::vector<double> velocities;

    for (const Particle& p : species.part_list) 
    {
        positions.push_back(p.pos[0]);
        positions.push_back(p.pos[1]);
        velocities.push_back(p.vel[0]);
        velocities.push_back(p.vel[1]);
    }

    dataset_pos.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_vel.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
}

//---------------------------

void Output::write_field_data(int ts)
{
    
    std::string pot_group_name  = "pot";
    std::string efieldx_group_name = "efieldx";
    std::string efieldy_group_name = "efieldy";

    Group pot_subgroup;
    Group efieldx_subgroup;
    Group efieldy_subgroup;

    if (!field_data_group.exists(pot_group_name))
    {
        // Subgroup doesn't exist, create it
        pot_subgroup = field_data_group.createGroup(pot_group_name);
        efieldx_subgroup = field_data_group.createGroup(efieldx_group_name);
        efieldy_subgroup = field_data_group.createGroup(efieldy_group_name);
    }
    else
    {
        // Subgroup already exists, retrieve it
        pot_subgroup = field_data_group.openGroup(pot_group_name);
        efieldx_subgroup = field_data_group.openGroup(efieldx_group_name);
        efieldy_subgroup = field_data_group.openGroup(efieldy_group_name);
    }
    
    //pot_subgroup = field_data_group.createGroup(subgroup_name);
    

    hsize_t nx = domain.nx;
    hsize_t ny = domain.ny;

    hsize_t dims_den[2] = {ny, nx};

    hsize_t Rank = 2;

    DataSpace dataspace_pot(Rank, dims_den);
    DataSpace dataspace_efieldx(Rank, dims_den);
    DataSpace dataspace_efieldy(Rank, dims_den);

    H5::DataSet dataset_pot = pot_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_pot);
    H5::DataSet dataset_efieldx = efieldx_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_efieldx);
    H5::DataSet dataset_efieldy = efieldy_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_efieldy);
    // Prepare data buffer

    std::vector<double> pot(nx * ny);
    std::vector<double> efieldx(nx * ny);
    std::vector<double> efieldy(nx * ny);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int j = 0; j < ny; ++j) // j is y-index (rows)
    {
        for (int i = 0; i < nx; ++i) // i is x-index (columns)
        {
            pot[j * nx + i] = domain.phi(i,j);
            efieldx[j * nx + i] = domain.efx(i,j);
            efieldy[j * nx + i] = domain.efy(i,j);
        }
    }

    // Write the density data to the dataset
    dataset_pot.write(pot.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_efieldx.write(efieldx.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_efieldy.write(efieldy.data(), H5::PredType::NATIVE_DOUBLE);
}

//-----------------den data----------------------------------------
void Output::write_den_data(int ts,  Species& species)
{
    
    std::string subgroup_name = "den_" + species.name;

    Group den_subgroup;

    // Check if the group already exists in the map
    auto it = den_subgroups.find(species.name);
    if (it != den_subgroups.end()) 
    {
        // Group already exists, retrieve it from the map
        den_subgroup = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        den_subgroup = field_data_group.createGroup(subgroup_name);
        den_subgroups[species.name] = den_subgroup;
    }
    
    hsize_t nx = domain.nx;
    hsize_t ny = domain.ny;

    hsize_t dims_den[2] = {ny, nx};

    hsize_t Rank = 2;

    DataSpace dataspace_den(Rank, dims_den);

    H5::DataSet dataset_den = den_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_den);

    // Prepare data buffer
    std::vector<double> density(nx * ny);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int j = 0; j < ny; ++j) 
    {
        for (int i = 0; i < nx; ++i)
        {
            density[j * nx + i] = species.den(i,j);
        }
    }

    // Write the density data to the dataset
    dataset_den.write(density.data(), H5::PredType::NATIVE_DOUBLE);
}



void Output::write_metadata(int nx, int ny, int NUM_TS, int write_int, int write_int_phase, double DT, double density, int save_fig, int normscheme, 
    int subcycleint, double LDe, double LDi, double wpe, double wpi,int spno, double GAS_DENSITY, double Bx, double By, double Bz)
{
    // Write metadata attributes within the metadata group
    //Group metadata_group = file.openGroup("/metadata");

    metadata_group.createAttribute("Nx", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &nx);
    metadata_group.createAttribute("Ny", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &ny);
    metadata_group.createAttribute("NUM_TS", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &NUM_TS);
    metadata_group.createAttribute("write_int", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int);
    metadata_group.createAttribute("write_int_phase", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int_phase);
    metadata_group.createAttribute("DT_coeff", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &DT);
    metadata_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &density);
    metadata_group.createAttribute("save_fig", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &save_fig);
    metadata_group.createAttribute("norm_scheme", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &normscheme);
    metadata_group.createAttribute("sub_cycle_interval", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &subcycleint);
    metadata_group.createAttribute("LDe", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &LDe);
    metadata_group.createAttribute("LDi", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &LDi);
    metadata_group.createAttribute("wpe", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &wpe);
    metadata_group.createAttribute("wpi", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &wpi);
    metadata_group.createAttribute("spno", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &spno);
    metadata_group.createAttribute("GAS_DENSITY", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &GAS_DENSITY);
    metadata_group.createAttribute("B", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &Bx);
    metadata_group.createAttribute("theta", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &By);
    metadata_group.createAttribute("azimuth", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &Bz);
    //metadata_group.createAttribute("max_ele_coll_freq", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &max_electron_collision_freq);
    metadata_group.close();
}

/*
void Output::write_species_metadata(std::vector<Species> &species_list)
{
    for (Species& sp : species_list)
    {
        std::string species_group_name = sp.name;
        Group species_group;
   
        species_group = metadata_species.createGroup(species_group_name);
    
        // Write attributes specific to the species
        species_group.createAttribute("name", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.name.c_str());
        species_group.createAttribute("mass", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.mass);
        species_group.createAttribute("charge", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.charge);
        species_group.createAttribute("spwt", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.spwt);
        species_group.createAttribute("temperature", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.temp);
        species_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.fract_den);
        species_group.createAttribute("num_particles", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &sp.numparticle);
        species_group.createAttribute("streaming_velocity", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.vsx);

        // Close the species group after writing the metadata
        species_group.close();
    }
}
*/

void Output::write_species_metadata(std::vector<Species>& species_list)
{
    // Create a vector to store species names in the correct order
    std::vector<std::string> species_order;

    for (Species& sp : species_list)
    {
        std::string species_group_name = sp.name;
        Group species_group;
   
        species_group = metadata_species.createGroup(species_group_name);
    
        // Write attributes specific to the species
        species_group.createAttribute("name", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.name.c_str());
        species_group.createAttribute("mass", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.mass);
        species_group.createAttribute("charge", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.charge);
        species_group.createAttribute("spwt", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.spwt);
        species_group.createAttribute("temperature", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.temp);
        species_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.fract_den);
        species_group.createAttribute("num_particles", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &sp.numparticle);
        species_group.createAttribute("streaming_velocity", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.vsx);

        // Add species name to species_order vector
        species_order.push_back(sp.name);

        // Close the species group after writing the metadata
        species_group.close();
    }

    // Create an attribute to store the species order as in species_list, use hsize_t array for the dimension
    hsize_t dims[1] = {species_order.size()};
    DataSpace order_space(H5S_SIMPLE, dims);

    // Define a variable-length string type for the species names as each species name have diffrent lenght strings
    StrType str_type(PredType::C_S1, H5T_VARIABLE);

    Attribute order_attr = metadata_species.createAttribute("species_order", str_type, order_space);

    std::vector<const char*> species_name_pointers;
    for (const auto& name : species_order)
    {
        species_name_pointers.push_back(name.c_str());
    }

    // Write the array of species names as a string array (variable-length strings)
    order_attr.write(str_type, species_name_pointers.data());
}

void Output::write_ke()
{
    // Define the name for the dataset
    std::string datasetName = "kinetic_energy";

    // Define the dimensions of the dataset
    hsize_t ny = int(domain.NUM_TS/domain.write_interval) + 1; //rows
    hsize_t nx = hsize_t(3*domain.species_no + 3); //column

    hsize_t dims_energy[2] = {ny, nx};

    // Create dataspace for the dataset
    DataSpace dataspace_energy(2, dims_energy);

    // Create the dataset within the time_group
    H5::DataSet dataset_energy = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace_energy);

    // Allocate memory for the data
    double* ke_data = new double[nx * ny];

    // Fill the data array
    for (hsize_t i = 0; i < ny; ++i)
    {
        for (hsize_t j = 0; j < nx; ++j)
        {
            ke_data[i * nx + j] = store_ke(i,j);
        }
    } 

    // Write the data to the dataset
    dataset_energy.write(ke_data, H5::PredType::NATIVE_DOUBLE);

    // Deallocate memory after writing
    //delete[] ke_data;
}
//*/


void Output::storeKE_to_matrix(int ts, std::vector<Species> &species_list)
{
    // Determine normalization species
    Species &norm_species = (domain.SolverType == "hybrid" || domain.SolverType == "gshybrid") ? species_list[0] : ((domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0]);

    // Timestep index for storage
    int k = ts / domain.write_interval;

    // Store timestep value
    store_ke(k, 0) = ts * domain.DT;

    //Store kinetic energy for each species (3 components each)
    int j = 1;  // Start after timestep column
    for (Species &sp : species_list)
    {
        vec<double> ke = sp.Compute_KE(norm_species);  // Returns KE vector [x, y, z]
        store_ke(k, j)     = ke(0);  // KE in x-direction
        store_ke(k, j + 1) = ke(1);  // KE in y-direction
        store_ke(k, j + 2) = ke(2);  // KE in z-direction
        j += 3;  // Increment by 3 for next species
    }

    // Store potential energy in the last two column
    vec<double> pe = domain.Compute_PE(norm_species);
    store_ke(k, j)     = pe(0);  // Potential energy Ex
    store_ke(k, j + 1) = pe(1);  // Potential energy Ey
}


void Output::diagnostics(int ts, std::vector<Species> &species_list, const PlotFlags &flags)
{
    
    if (domain.diagtype == "off")
    {
        std::cout << "TS: " << ts << std::endl;
        return; // Exit the function if diagnostics are turned off
    }

    Species& norm_species = (domain.SolverType == "hybrid" || domain.SolverType == "gshybrid") ? species_list[0]:((domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0]);

    if (domain.diagtype == "basic")
    {

        vec<double> phivec = flatten(domain.phi);
        double max_phi = phivec(0);

        for (int i = 0; i < int(domain.nx * domain.ny); i++)
        {
            if (phivec(i) > max_phi)
            {
                max_phi = phivec(i);
            }
        }
        
        std::cout << "TS: " << ts << "\t";

        if(domain.SolverType != "spectral")
        {
            std::cout <<" norm: " <<std::scientific<<domain.error;
        }
        std::cout <<" delta_phi: " << std::fixed << std::setprecision(2) << (max_phi - phivec(0));

        if (domain.bc == "open")
        {
            for (Species& sp : species_list)
            {
                std::cout << " n_" << std::setw(4) << sp.name << ":" << sp.part_list.size();
            }   
        }

        std::cout<<endl;
    }

    if (domain.diagtype == "full")
    {
        
        if (domain.normscheme == 5)
        {
            std::cout << std::scientific << std::setprecision(precision);
        }
        else
        {
            std::cout << std::fixed << std::setprecision(precision);
        }
        
        vec<double> phivec = flatten(domain.phi);
        double max_phi = phivec(0);

        for (int i = 0; i < int(domain.nx * domain.ny); i++)
        {
            if (phivec(i) > max_phi)
            {
                max_phi = phivec(i);
            }
        }
        
        std::cout << "TS: " << ts << "\t";

        if(domain.SolverType != "spectral")
        {
            std::cout <<" norm: " <<std::scientific<<domain.error;
        }
        std::cout <<" delta_phi: "  << (max_phi - phivec(0));

        if (domain.bc == "open")
        {
            for (Species& sp : species_list)
            {
                std::cout << " n_" << std::setw(4) << sp.name << ":" << sp.part_list.size();
            }   
        }
    
        double total_kinetic_energy = 0.0;
        double potential_energy_value = 0.0;
        vec<double> total_ke_components(3);

        if (flags.ke_components == 1 || flags.total_energy == 1)
        {
            for (Species& sp : species_list)
            {
                vec<double> ke = sp.Compute_KE(norm_species);
                total_ke_components += ke;
                total_kinetic_energy += ke(0) + ke(1) + ke(2);
                if(domain.bc == "pbc")
                {
                    std::cout << " KE_" << std::setw(4) << sp.name << ": " << ke(0) + ke(1) + ke(2);
                }
                
            }

            //now there are two pe component efx**2 and efy**2
            vec<double> pe_comp = domain.Compute_PE(norm_species);
            potential_energy_value = pe_comp(0) + pe_comp(1);
            
            if(flags.total_energy == 1)
            {
                std::cout << " PEx: " << pe_comp(0);
                std::cout << " PEy: " << pe_comp(1);
                std::cout << " PE : "  << potential_energy_value;
                std::cout << " KE : "  << total_kinetic_energy;
                std::cout << " TE : "  << total_kinetic_energy + potential_energy_value;
            }
            
            vec<double> total_momentum(3);
            total_momentum = 0;
            for (Species& sp : species_list)
            {
                total_momentum += sp.Compute_Momentum(norm_species);
            }
            
            /*
            if(domain.bc == "pbc")
            {
                std::cout << " p_x: " << total_momentum(0);
                std::cout << " p_y: " << total_momentum(1);
                std::cout << " p_z: " << total_momentum(2);
                //std::cout << " p: " << std::fixed << std::setprecision(precision) 
            //<< sqrt(total_momentum(0)*total_momentum(0) + total_momentum(1)*total_momentum(1) + total_momentum(2)*total_momentum(2));
            }*/
            
            // Energy history store
            time_steps.push_back(static_cast<double>(ts));
            kinetic_energy.push_back(total_kinetic_energy);
            potential_energy.push_back(potential_energy_value);
            total_energy.push_back(total_kinetic_energy + potential_energy_value);
            Ke_x.push_back(total_ke_components(0));
            Ke_y.push_back(total_ke_components(1));
            Ke_z.push_back(total_ke_components(2));
        }

        std::cout << std::endl;

        plt::ion();

        // Phase Space Plot (x vs vx)
        if ((flags.phase_space == 1 || flags.phase_space == 2) && flags.species_index < species_list.size())
        {
            std::vector<double> x, vx;
            int N = species_list[flags.species_index].part_list.size();
            for (int i = 0; i < N; ++i)
            {
                (flags.phase_space == 1) ? x.push_back(species_list[flags.species_index].part_list[i].pos[0]):
                x.push_back(species_list[flags.species_index].part_list[i].pos[1]);
                (flags.phase_space == 1) ? vx.push_back(species_list[flags.species_index].part_list[i].vel[0]):
                vx.push_back(species_list[flags.species_index].part_list[i].vel[1]);
            }

            plt::figure(1);
            plt::clf();
            plt::scatter(x, vx, 1.0);
            (flags.phase_space == 1) ? plt::xlabel("x"): plt::xlabel("y");
            (flags.phase_space == 1) ? plt::ylabel("vx"): plt::ylabel("vy");
            //plt::xlim(0, domain.nx);
            //plt::ylim(0, domain.ny);
        }

        // Config Space Plot (x vs y)
        if (flags.config_space == 1 && flags.species_index < species_list.size())
        {
            std::vector<double> x, y;
            int N = species_list[flags.species_index].part_list.size();
            for (int i = 0; i < N; ++i)
            {
                x.push_back(species_list[flags.species_index].part_list[i].pos[0]);
                y.push_back(species_list[flags.species_index].part_list[i].pos[1]);
            }

            plt::figure(2);
            plt::clf();
            plt::scatter(x, y, 1.0);
            plt::xlabel("x");
            plt::ylabel("y");
            plt::xlim(0, domain.nx);
            plt::ylim(0, domain.ny);
        }


        // Electric Field Plot
        if (flags.electric_field == 1 || flags.electric_field == 2)
        {
            // Convert domain.efx to std::vector<std::vector<double>>
            //std::vector<std::vector<double>> phi_data(domain.nx, std::vector<double>(domain.ny));
            std::vector<std::vector<double>> efield_data(domain.ny, std::vector<double>(domain.nx));

            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    efield_data[j][i] = (flags.electric_field == 1) ? domain.efx(i, j) : domain.efy(i, j);
                }
            }

            // Plot using imshow
            plt::figure(3);
            plt::clf();
            plt::imshow(efield_data, "coolwarm", "lower","bilinear");
            //plt::colorbar();
            (flags.electric_field == 1) ? plt::title("Electric Field E_x") : plt::title("Electric Field E_y");
            //plt::title("Electric Field E");
            plt::xlabel("x");
            plt::ylabel("y");
        }
        
        // KE Components Plot
        if (flags.ke_components == 1)
        {
            plt::figure(4);
            plt::clf();
            plt::plot(time_steps, Ke_x, {{"label", "KE_x"}, {"color", "blue"}});
            plt::plot(time_steps, Ke_y, {{"label", "KE_y"}, {"color", "red"}});
            plt::plot(time_steps, Ke_z, {{"label", "KE_z"}, {"color", "green"}});
            plt::title("Kinetic Energy Components vs Time");
            plt::xlabel("Time Step");
            plt::ylabel("Kinetic Energy component");
            plt::legend({{"loc", "upper right"}});
        }

        // Total Energy Plot
        if (flags.total_energy == 1)
        {
            plt::figure(5);
            plt::clf();
            plt::plot(time_steps, kinetic_energy, {{"label", "Total Kinetic Energy"}, {"color", "blue"}});
            plt::plot(time_steps, potential_energy, {{"label", "Potential Energy"}, {"color", "red"}});
            plt::plot(time_steps, total_energy, {{"label", "Total Energy"}, {"color", "green"}});
            plt::title("Total Energy vs Time");
            plt::xlabel("Time Step");
            plt::ylabel("Energy");
            plt::legend({{"loc", "upper right"}});
        }

        
        // Potential Field Plot (phi) using imshow
        if (flags.potential_field == 1)
        {
            // Convert domain.phi to std::vector<std::vector<double>>
            //std::vector<std::vector<double>> phi_data(domain.nx, std::vector<double>(domain.ny));
            std::vector<std::vector<double>> phi_data(domain.ny, std::vector<double>(domain.nx));

            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    phi_data[j][i] = domain.phi(i, j);
                }
            }

            // Plot using imshow
            plt::figure(6);
            plt::clf();
            plt::imshow(phi_data, "coolwarm", "lower","bilinear");
            //plt::colorbar();
            plt::title("Potential Field φ");
            plt::xlabel("x");
            plt::ylabel("y");


            // Prepare Dirichlet boundary points for scatter plot
            std::vector<double> x_points;
            std::vector<double> y_points;

            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    if (domain.dirichlet_nodes(i, j))  // Assuming this is your boolean matrix
                    {
                        x_points.push_back(i);
                        y_points.push_back(j);
                    }
                }
            }

            // Plot Dirichlet BC as black dots
            plt::scatter(x_points, y_points, 5.0, {{"color", "black"}});

            plt::title("potential(with dirichlet nodes)");
            plt::xlabel("x");
            plt::ylabel("y");
        }

        // density contour plot
        if (flags.density_contour == 1)
        {
            // Convert domain.phi to std::vector<std::vector<double>>
            //std::vector<std::vector<double>> phi_data(domain.nx, std::vector<double>(domain.ny));
            std::vector<std::vector<double>> density_data(domain.ny, std::vector<double>(domain.nx));

            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    density_data[j][i] = species_list[flags.species_index].den(i,j);
                }
            }

            // Plot using imshow
            plt::figure(7);
            plt::clf();
            plt::imshow(density_data, "coolwarm", "lower","bilinear");
            //plt::colorbar();
            plt::title("Density");
            plt::xlabel("x");
            plt::ylabel("y");
        }


        plt::pause(0.1);
        plt::show();

    }
 
}
