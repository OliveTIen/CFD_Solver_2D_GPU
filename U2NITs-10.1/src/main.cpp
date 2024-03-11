#include "../include/main.hpp"

int main(int argc, char* argv[]) {
	const int MAX_STRING_SIZE = 1024;
	int nZone = 1;

	char config_file_name[MAX_STRING_SIZE];
	bool dry_run = false;
	int num_threads = omp_get_max_threads();
	bool use_thread_mult = false;
	std::string filename = "default.cfg";

	CLI::App app{"SU2 v8.0.0 \"Harrier\", The Open-Source CFD Code"};
	app.add_flag("-d,--dryrun", dry_run, "Enable dry run mode.\n"
		"Only execute preprocessing steps using a dummy geometry.");
	app.add_option("-t,--threads", num_threads, "Number of OpenMP threads per MPI rank.");
	app.add_flag("--thread_multiple", use_thread_mult, "Request MPI_THREAD_MULTIPLE thread support.");
	app.add_option("configfile", filename, "A config file.")->check(CLI::ExistingFile);

	CLI11_PARSE(app, argc, argv)

	CDriver* driver = nullptr;

	if (dry_run) {

		/*--- Dry Run. ---*/
		driver = new CDummyDriver(config_file_name, nZone);

	}

	driver->StartSolver();

	/*--- Finalize solver, delete all the containers, close history file, exit SU2. ---*/

	driver->Finalize();

	delete driver;

	return EXIT_SUCCESS;
}

/*
PS D:\tgl\Local\HPC\U2NITS-10.0\x64\Debug> ./U2NITs-10.1.exe -h
SU2 v8.0.0 "Harrier", The Open-Source CFD Code
Usage: D:\tgl\Local\HPC\U2NITS-10.0\x64\Debug\U2NITs-10.1.exe [OPTIONS] [configfile]

Positionals:
  configfile TEXT:FILE        A config file.

Options:
  -h,--help                   Print this help message and exit
  -d,--dryrun                 Enable dry run mode.
							  Only execute preprocessing steps using a dummy geometry.
  -t,--threads INT            Number of OpenMP threads per MPI rank.
  --thread_multiple           Request MPI_THREAD_MULTIPLE thread support.


*/