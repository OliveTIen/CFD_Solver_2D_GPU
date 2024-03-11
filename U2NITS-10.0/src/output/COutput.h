
#ifndef _COUTPUT_H_
#define _COUTPUT_H_

#include "ConsolePrinter.h"
#include "FieldWriter.h"
#include "HistWriter.h"
#include "LogWriter.h"
#include "ResidualCalculator.h"
namespace U2NITS {
	class COutput {
	private:
		//ConsolePrinter* m_pConsolePrinter;
		//FieldWriter* m_pFieldWriter;
		HistWriter* m_pHistWriter;
		//LogWriter* m_pLogWriter;
		//ResidualCalculator* m_pResidualCalculator;
	public:
		COutput(std::string histFilePath) {
			//m_pConsolePrinter = new ConsolePrinter();
			//m_pFieldWriter = new FieldWriter();
			m_pHistWriter = new HistWriter(histFilePath);
			//m_pLogWriter = new LogWriter();
			//m_pResidualCalculator = new ResidualCalculator();
		}
		~COutput() {
			delete m_pHistWriter;
		}


	};
}

#endif // _COUTPUT_H_