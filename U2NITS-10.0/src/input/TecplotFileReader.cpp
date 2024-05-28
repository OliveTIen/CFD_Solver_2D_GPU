#include "TecplotFileReader.h"

TecplotFileReader* TecplotFileReader::getInstance() {

    if (p_instance == nullptr) {
        p_instance = new TecplotFileReader();
    }
    return p_instance;
}

TecplotFileReader* TecplotFileReader::p_instance = nullptr;
