/*
�����������Ե��࣬�ʺ����Ž�ģʽ
*/

#ifndef _CSOLVER_H_
#define _CSOLVER_H_
#include "equation/CEquation.h"
#include "platform/CPlatform.h"

namespace U2NITS {
	class CSolver {
	private:
		//CEquation* m_equation;
		//CPlatform* m_platform;
	public:
		virtual void initialize() = 0;
		virtual void iterate() = 0;
		virtual void updateResidual() = 0;
		virtual void finalize() = 0;
	};
}
#endif // !_CSOLVER_H_

/*
���裺���Ż���
��ȡ���ã����Ʋ����� FileReader.readConfig()
��ʼ����������ȡ�����ļ� / ��ȡ�����ļ� + ��ʼ�� FileReader.readField()
�����ʼ��Ϣ��������Ϣ���߽�����ȣ�ConsolePrinter.updateScreen()
����в��ļ�ͷ FileIOManager.writeHistFileHead()
��ʱ����ʼ
��ʼ����Դ��GPU�ڴ棩Solver.initialize()
������ʼ
���²������ݣ��ļ�����
���㣨������ Solver.iterate()
���㣨�в�������ʱ���� Solver.updateResidual()
���������Ϣ��������¼��һ�ι��λ�á�������ȣ�ConsolePrinter.updateScreen()
����ļ����������вFileWriter.write()
�ж��Ƿ�����ѭ������ֹ��this->checkStop()
���������Ϣ����ʾ�ʣ�ConsolePrinter.updateScreen()
�ͷ���Դ Solver.finalize()
*/
