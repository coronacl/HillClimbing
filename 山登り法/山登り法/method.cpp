#include "stdafx.h"
#include "gnuplot.h"

using namespace gnuplot;

//System::String��std::string�֕ϊ�
void mw::String2string(System::String ^ s, std::string& os) {
	using namespace Runtime::InteropServices;
	const char* chars =
		(const char*)(Marshal::StringToHGlobalAnsi(s)).ToPointer();
	os = chars;
	Marshal::FreeHGlobal(IntPtr((void*)chars));
}

//System::String^ mw::string2String(std::string str)
//{
//	System::String^ dst = gcnew System::String(str.c_str());
//	return(dst);
//}

//int��System::String�֕ϊ�
System::String^ mw::int2String(int tmpi)
{
	std::string str = to_string(tmpi);
	System::String^ dst = gcnew System::String(str.c_str());
	return(dst);
}

//double��System::String�֕ϊ�
System::String^ mw::double2String(double tmpd)
{
	std::string str = to_string(tmpd);
	System::String^ dst = gcnew System::String(str.c_str());
	return(dst);
}
//�t�@�C���܂ł̐�΃p�X���擾
System::String^ mw::getFilePath()
{
	////�r�b�g�}�b�v�t�@�C���̎擾 
	//̧�ق��J���_�C�A���O�̍쐬 
	OpenFileDialog^ dlg = gcnew OpenFileDialog;

	//�t�@�C���t�B���^ 
	//dlg->Filter = "�摜̧��(*.bmp,*.jpg,*.png,*.tif,jp2)|*.bmp;*.jpg;*.png;*.tif;*.jp2";

	//�_�C�A���O�̕\�� 
	if (dlg->ShowDialog() == System::Windows::Forms::DialogResult::Cancel){
		cout << "cannot open file" << endl;
	}

	//System::String^�^�̃t�@�C����
	System::String^ strFilename = dlg->FileName;
	return(strFilename);
}

//���̃t�@�C���܂ł̃p�X���擾(�t�@�C�����͊܂܂Ȃ�)
System::String^ mw::getFileListPath()
{
	////�r�b�g�}�b�v�t�@�C���̎擾 
	//̧�ق��J���_�C�A���O�̍쐬 
	OpenFileDialog^ dlg = gcnew OpenFileDialog;

	//�t�@�C���t�B���^ 
	dlg->Filter = "dcm files (*.dcm)|*.dcm";

	//�_�C�A���O�̕\�� 
	if (dlg->ShowDialog() == System::Windows::Forms::DialogResult::Cancel){
		MessageBox::Show("dicom�t�@�C����I��ł�");
	}

	//System::String^�^�̃t�@�C����
	System::String^ strFilePath = dlg->FileName;//�t�@�C�����܂ł̐�΃p�X�擾
	System::String^ strFileName = dlg->SafeFileName;//�t�@�C�����擾
	int num_name = strFileName->Length;
	int num_full = strFilePath->Length;
	int num_dis = num_full - num_name;
	System::String^ dst = strFilePath->Remove(num_dis);

	return(dst);
}

//���̃f�B���N�g���̃t�@�C�����ׂĎ擾
array<System::String^>^ mw::getFileList()
{
	System::String^ strFilename = getFileListPath();
	array<System::String^>^ dir = Directory::GetFiles(strFilename);

	//cout << "name of files = " << endl;;
	for (int i = 0; i < dir->Length; i++){
		//Console::WriteLine(dir[i]);
		//cout << endl;
		if (!dir[i]->Contains(".dcm")){
			MessageBox::Show("�t�H���_����dcm�݂̂ɂ��Ă�������");
		}

	}

	return(dir);
}

//Sytem::String��std::string�֕ϊ�
void mw::MarshalString(System::String ^ s, string& os) {
	using namespace Runtime::InteropServices;
	const char* chars =
		(const char*)(Marshal::StringToHGlobalAnsi(s)).ToPointer();
	os = chars;
	Marshal::FreeHGlobal(IntPtr((void*)chars));
}

//int��std::string�֕ϊ�
string mw::int2string(int number)
{
	return(to_string(number));
}




////csv�t�@�C��������
//(template�֐�)
//template<typename T_n>
//void makeCsv(cv::Mat src)
//{
//	FILE *fp;
//	char filename[] = "dataByEachBIN3gauss.csv";
//
//	errno_t error;
//	if (error = fopen_s(&fp, filename, "w") != 0){
//		fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
//	}
//
//	T_n tmp;
//
//	for (int j = 0; j < src.rows; ++j){
//		for (int i = 0; i < src.cols; ++i){
//			//tmp = PIXd(src, i, j);
//			tmp = src.at<T_n>(j, i);
//			fprintf(fp, "%lf", tmp);
//			if (i == src.cols - 1){
//				fprintf(fp, "\n", tmp);
//			}
//			else{
//				fprintf(fp, ",");
//			}
//		}
//	}
//	fclose(fp);
//}

Mat_<float> mw::csv2matF(string file_name)
{
	cout << "loading...";
	//�s���Ɨ񐔂̃`�F�b�N--------------------------------
	ifstream tmp_ifs(file_name);
	if (tmp_ifs.fail()){
		cout << "I cannot open " << file_name << endl;
	}

	int num_rows = 0, num_cols = 0;
	string tmp_str, tmp_segment;
	bool flag = true;
	while (getline(tmp_ifs, tmp_str)){
		num_rows++;
		if (flag){
			std::istringstream tmp_stream(tmp_str);
			while (getline(tmp_stream, tmp_segment, ',')){
				num_cols++;
			}
			flag = false;
		}
	}
	//cout << "rows = " << num_rows << " cols = " << num_cols << endl;
	tmp_ifs.close();

	//Mat�ւ̓���---------------------------------------
	static Mat dst = Mat::zeros(num_rows, num_cols, CV_32FC1);
	ifstream ifs(file_name);
	string str;
	string seg;
	const char* tmpc;
	int i = 0, j = 0;
	while (getline(ifs, str))
	{
		std::istringstream stream(str);
		while (getline(stream, seg, ','))
		{
			//cout <<"(" << i << "," << j << ") " << seg.c_str() << endl;
			tmpc = seg.c_str();
			dst.at<float>(j, i) = atof(tmpc);
			i++;
		}
		i = 0;
		j++;
	}
	//cout << dst << endl;
	cout << "end" << endl;
	return(dst);
}

//�^��uchar�݂̂���
Mat& mw::csv2mat(string file_name)
{
	cout << "loading...";
	//�s���Ɨ񐔂̃`�F�b�N--------------------------------
	ifstream tmp_ifs(file_name);
	if (tmp_ifs.fail()){
		cout << "I cannot open " << file_name << endl;
	}

	int num_rows = 0, num_cols = 0;
	string tmp_str, tmp_segment;
	bool flag = true;
	while (getline(tmp_ifs, tmp_str)){
		num_rows++;
		if (flag){
			std::istringstream tmp_stream(tmp_str);
			while (getline(tmp_stream, tmp_segment, ',')){
				num_cols++;
			}
			flag = false;
		}
	}
	//cout << "rows = " << num_rows << " cols = " << num_cols << endl;
	tmp_ifs.close();

	//Mat�ւ̓���---------------------------------------
	static Mat dst = Mat::zeros(num_rows, num_cols, CV_8UC1);
	ifstream ifs(file_name);
	string str;
	string seg;
	const char* tmpc;
	int i = 0, j = 0;
	while (getline(ifs, str))
	{
		std::istringstream stream(str);
		while (getline(stream, seg, ','))
		{
			//cout <<"(" << i << "," << j << ") " << seg.c_str() << endl;
			tmpc = seg.c_str();
			dst.at<uchar>(j, i) = atoi(tmpc);
			i++;
		}
		i = 0;
		j++;
	}
	//cout << dst << endl;
	cout << "end" << endl;
	return(dst);
}




//�񎟌��O���t������
void mw::makeGraph2d(cv::Mat hist)
{
	vector<double> x, y;

	for (int i = 0; i < hist.cols; ++i){
		x.push_back(i);
	};

	double tmpd;
	for (int j = 0; j < hist.cols; ++j){
		tmpd = (double)hist.at<int>(0, j);
		y.push_back(tmpd);
	};


	static CGnuplot gp;
	gp.SetLabel("x", "y");
	//gp.SetView(33, 71, 1, 1);//���_
	gp.Plot(x, y);//�v���b�g�͍Ō�
}

//�O�����O���t��`��
void mw::makeGraph3d(cv::Mat src)
{
	vector<double> x, y;
	vector<vector<double> > z;

	for (int i = 0; i < src.cols; ++i){ x.push_back(i); };
	for (int j = 0; j < src.rows; ++j){ y.push_back(j); };

	double tmp;
	for (int j = 0; j < src.rows; ++j){
		vector<double> _z;
		for (int i = 0; i < src.cols; ++i){
			tmp = (double)src.at<double>(Point(i, j));
			_z.push_back(tmp);
		}
		z.push_back(_z);
	}

	static CGnuplot* gp = new CGnuplot();
	gp->SetLabel("y", "x");
	gp->SetView(0.0, 90.0, 1.0, 1.0);
	gp->Plot(x, y, z);

	//__KEYWAIT__;

	//gp->~CGnuplot();

}

void mw::swqp(int& x, int& y)
{
	int tmp;
	tmp = x;
	x = y;
	y = tmp;
}

void mw::makeRandom(int* ans)
{
	double tmp;
	int index;
	for (int i = 0; i < 256; i++){
		tmp = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
		index = (int)(tmp * 256);
		swap(ans[i], ans[index]);
		//cout << i << "::" << index << endl;
	}
}

Point mw::calGlavityCenter(Mat_<double> src)
{
	double sum_dis = 0.0;
	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){
			sum_dis +=  src.at<double>(j,i);
		}
	}

	double disX = 0.0, disY = 0.0;
	double dis;
	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){
			dis =  src.at<double>(j, i) ;
			disX += dis / sum_dis * i;
			disY += dis / sum_dis * j;
		}
	}
	return(Point((int)(disX + 0.5), (int)(disY+ 0.5)));
}

//ZNCC���v�Z�����
double mw::calZNCC(Mat_<int>& bf, Mat_<int>& af)
{
	if (bf.rows != af.rows || bf.cols != af.cols){
		std::cout << "�摜�T�C�Y���Ⴂ�܂�" << endl;
		return 0;
	}

	int sum_bf = 0, sum_af = 0;
	for (int j = 0; j < bf.rows; ++j){
		for (int i = 0; i < bf.cols; ++i){
			sum_bf = sum_bf + PIXi(bf, i, j);
			sum_af = sum_af + PIXi(af, i, j);
		}
	}
	double ave_bf = (double)sum_bf / (bf.rows * bf.cols);
	double ave_af = (double)sum_af / (af.rows * af.cols);

	double zncc = 0.0;
	double numerator = 0.0;//���q
	double denominator = 0.0;//����
	double tmpA = 0.0;
	double tmpB = 0.0;
	for (int j = 0; j < bf.rows; ++j){
		for (int i = 0; i < bf.cols; ++i){
			numerator += (PIXi(bf, i, j) - ave_bf) *(PIXi(af, i, j) - ave_af);
			tmpA = tmpA + (PIXi(bf, i, j) - ave_bf) * (PIXi(bf, i, j) - ave_bf);
			tmpB = tmpB + (PIXi(af, i, j) - ave_af) * (PIXi(af, i, j) - ave_af);
		}
	}
	denominator = sqrt(tmpA * tmpB);

	zncc = numerator / denominator;

	return (zncc);
}

////NMI�̌v�Z
//���S�ɏd�݂����A�q�X�g�O�������v�Z
double mw::calNMIinMask(Mat_<int> bf, Mat_<int> af, Mat_<int> mask, int BIN_DOWN)
{

	if (bf.rows != af.rows || bf.cols != af.cols || mask.rows != bf.cols){
		std::cout << "�摜�T�C�Y���Ⴂ�܂�" << endl;
		return 0;
	}

	//�؂������̈�̊K����ύX
	double max_val, min_val;
	Mat_<int> picNMI = Mat_<int>::zeros(bf.rows, bf.cols);
	cv::minMaxLoc(bf, &min_val, &max_val, NULL, NULL);
	for (int j = 0; j < bf.rows; j++){
		for (int i = 0; i < bf.cols; i++){
			picNMI[j][i] = (int)(((bf.at<int>(Point(i, j)) - min_val) / (max_val - min_val)) * (BIN_DOWN - 1));
		}
	}

	Mat_<int> picNMI_af = Mat_<int>::zeros(af.rows, af.cols);
	cv::minMaxLoc(af, &min_val, &max_val, NULL, NULL);
	for (int j = 0; j < af.rows; j++){
		for (int i = 0; i < af.cols; i++){
			picNMI_af[j][i] = (int)(((af.at<int>(Point(i, j)) - min_val) / (max_val - min_val)) * (BIN_DOWN - 1));
		}
	}

	int* hg_bf = new int[BIN_DOWN];
	int* hg_af = new int[BIN_DOWN];
	int** hg_2d;
	hg_2d = new int*[BIN_DOWN];
	for (int count = 0; count < BIN_DOWN; count++){
		hg_2d[count] = new int[BIN_DOWN];
	}

	for (int x = 0; x < BIN_DOWN; x++){
		hg_bf[x] = 0;
		hg_af[x] = 0;
	}
	for (int y = 0; y < BIN_DOWN; y++){
		for (int x = 0; x < BIN_DOWN; x++){
			hg_2d[y][x] = 0;
		}
	}

	int posX, posY;
	for (int y = 0; y < af.rows; y++){
		for (int x = 0; x < af.cols; x++){
			//posX = PIXi(bf, x, y) ;
			posX = picNMI.at<int>(y, x);
			hg_bf[posX] += 1 * mask.at<int>(y,x);

			//posY = PIXi(af, x, y) ;
			posY = picNMI_af.at<int>(y, x);
			hg_af[posY] += 1 * mask.at<int>(y, x);

			hg_2d[posY][posX] += 1 * mask.at<int>(y, x);

		}
	}
	////�K�E�V�A���t�B���^�̓K�p�o�[�W����
	//
	Mat_<double> hg_bf2 = Mat_<double>::zeros(BIN_DOWN, 1);
	Mat_<double> hg_af2 = Mat_<double>::zeros(BIN_DOWN, 1);
	Mat_<double> hg_2d2 = Mat_<double>::zeros(BIN_DOWN, BIN_DOWN);
	for (int i = 0; i < BIN_DOWN; i++){
		hg_bf2.at<double>(i, 0) = hg_bf[i];
		hg_af2.at<double>(i, 0) = hg_af[i];
	}
	for (int y = 0; y < BIN_DOWN; y++){
		for (int x = 0; x < BIN_DOWN; x++){
			hg_2d2.at<double>(y, x) = hg_2d[y][x];
		}
	}
	int gauss_size = 5;
	int gauss_sigma = 3;
	//cv::GaussianBlur(hg_bf2, hg_bf2, Size(gauss_size, gauss_size), gauss_sigma, gauss_sigma);
	//cv::GaussianBlur(hg_af2, hg_af2, Size(gauss_size, gauss_size), gauss_sigma, gauss_sigma);//�����Ⴄ�C�����邯�ǂƂ肠�������u
	//cv::GaussianBlur(hg_2d2, hg_2d2, Size(gauss_size, gauss_size), gauss_sigma, gauss_sigma);

	//�ȉ��A�q�X�g�O��������NMI�̌v�Z
	//
	const double SUM_PIX = af.cols * af.rows;//��f��

	double* P_bf = new double[BIN_DOWN];
	double* P_af = new double[BIN_DOWN];
	double** P_2d;
	P_2d = new double*[BIN_DOWN];
	for (int count = 0; count < BIN_DOWN; count++){
		P_2d[count] = new double[BIN_DOWN];
	}

	for (int x = 0; x < BIN_DOWN; x++){
		P_bf[x] = hg_bf2.at<double>(x, 0) / SUM_PIX;
		P_af[x] = hg_af2.at<double>(x, 0) / SUM_PIX;
	}
	for (int y = 0; y < BIN_DOWN; y++){
		for (int x = 0; x < BIN_DOWN; x++){
			P_2d[y][x] = hg_2d2.at<double>(y, x) / SUM_PIX;
		}
	}

	double bfI = 0.0, afI = 0.0, nmi = 0.0;
	for (int x = 0; x < BIN_DOWN; x++){
		if (P_bf[x] != 0){
			bfI += -P_bf[x] * log(P_bf[x]);
		}
		if (P_af[x] != 0){
			afI += -P_af[x] * log(P_af[x]);
		}
	}

	double tmpX, tmpY, tmpXY;
	for (int y = 0; y < BIN_DOWN; y++){
		tmpY = P_af[y];
		for (int x = 0; x < BIN_DOWN; x++){
			tmpX = P_bf[x];
			tmpXY = P_2d[y][x];
			if (tmpX != 0 && tmpY != 0 && tmpXY != 0){
				nmi += tmpXY * log(tmpXY / tmpX / tmpY);
			}
		}
	}
	nmi = nmi / (bfI + afI) * 2.0;

	//�ȉ��������
	delete[] hg_bf;
	delete[] hg_af;
	delete[] P_bf;
	delete[] P_af;
	for (int i = 0; i<BIN_DOWN; i++) {
		delete[] hg_2d[i];
	}
	delete[] hg_2d;
	for (int i = 0; i<BIN_DOWN; i++) {
		delete[] P_2d[i];
	}
	delete[] P_2d;

	return(nmi);
}

//�s�[�N���ӂƂ��̑���SN��
double mw::calSN(Mat_<double> src)
{
	Mat ave;
	reduce(src, ave, 0, CV_REDUCE_AVG);//�s�ɕϊ�
	reduce(ave, ave, 1, CV_REDUCE_AVG);//��ɕϊ�
	//cout << "ave = " << ave << endl;
	double average = ave.at<double>(0, 0);

	double max_val, min_val;
	cv::Point2i max_loc;
	cv::minMaxLoc(src, &min_val, &max_val, NULL, &max_loc);

	//cout << "max=" << max_val << "max's location = " << max_loc << endl;

	const int SUM_PIX = src.rows * src.cols;
	double signal_sum = 0.0, noise_sum = 0.0;
	int num_signal = 0, num_noise = 0;
	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){
			if (src.at<double>(j, i) > average){
				noise_sum += src.at<double>(j, i);
				num_noise++;
			}
		}
	}

	const int mountain_size = 1;
	for (int j = max_loc.y - mountain_size; j <= max_loc.y + mountain_size; j++){
		for (int i = max_loc.x - mountain_size; i <= max_loc.x + mountain_size; i++){
			if (0 <= i && i < src.cols && 0 <= j && j < src.rows){
				if (src.at<double>(j, i) > average){
					noise_sum -= src.at<double>(j, i);
					signal_sum += src.at<double>(j, i);
					num_noise--;
					num_signal++;
				}
			}
		}
	}

	double signal, noise, sn;
	signal = (double)signal_sum / num_signal;
	noise = (double)noise_sum / num_noise;
	sn = signal / noise;

	//cout << "noise num = " << num_noise << endl;
	//cout << "signal num " << "= " << num_signal << endl;
	//cout << "signal=" << signal << " noise=" << noise << endl;
	//std::cout << "sn=" << sn << endl;

	return(sn);
}

//�ő�s�[�N�Ǝ��̃s�[�N��SN��
//�������A�R�̃T�C�Y��kernel_size*2+1�Ƃ��Ă���
double mw::calSN_ByPeak(Mat_<double>& src)
{
	double max_val;
	Point max_loc;
	const int kernel_size = 2;//���ۂ̃T�C�Y��kernel_size*2+1�Ƃ��A�����Ƃ���
	Mat_<double> kernel = Mat_<double>::zeros(kernel_size * 2 + 1, kernel_size * 2 + 1);
	Mat_<double> peak_table = Mat_<double>::zeros(src.rows, src.cols);

	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){

			if (0 <= i - kernel_size && i + kernel_size < src.cols){
				if (0 <= j - kernel_size && j + kernel_size < src.rows){
					kernel = src(Rect(i - kernel_size, j - kernel_size, kernel_size * 2 + 1, kernel_size * 2 + 1));
					cv::minMaxLoc(kernel, NULL, &max_val, NULL, &max_loc);
					if (max_loc.x == kernel_size && max_loc.y == kernel_size){
						peak_table.at<double>(j, i) = max_val;
						//cout << "find max = " << max_val << endl;
					}
				}
			}

		}
	}
	double first_peak = 0.0;
	double second_peak = 0.0;
	Point max_peak_loc;
	minMaxLoc(peak_table, NULL, &first_peak, NULL, &max_peak_loc);
	//cout << "loc.x = " << max_peak_loc.x << ",loc.y = " << max_peak_loc.y << endl;
	peak_table.at<double>(max_peak_loc.y, max_peak_loc.x) = 0.0;
	minMaxLoc(peak_table, NULL, &second_peak, NULL, &max_peak_loc);
	//cout << "loc.x = " << max_peak_loc.x << ",loc.y = " << max_peak_loc.y << endl;

	//�O��������O
	if (first_peak == 0.0){
		return(0.0);
	}
	if (first_peak != 0.0 && second_peak == 0.0){
		return(-100.0);
	}

	double sn = first_peak / second_peak;
	//cout << "sn = " << sn << endl;

	return(sn);
}

//�ő�s�[�N�Ǝ��̃s�[�N��SN��
//�������A�R�̃T�C�Y��kernel_size*2+1�Ƃ��Ă���
double mw::calSN_ByPeak(Mat_<double>& src, int mountain_r)
{
	double max_val;
	Point max_loc;
	int kernel_size = mountain_r;//���ۂ̃T�C�Y��kernel_size*2+1�Ƃ��A�����Ƃ���
	Mat_<double> kernel = Mat_<double>::zeros(kernel_size * 2 + 1, kernel_size * 2 + 1);
	Mat_<double> peak_table = Mat_<double>::zeros(src.rows, src.cols);

	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){

			if (0 <= i - kernel_size && i + kernel_size < src.cols){
				if (0 <= j - kernel_size && j + kernel_size < src.rows){
					kernel = src(Rect(i - kernel_size, j - kernel_size, kernel_size * 2 + 1, kernel_size * 2 + 1));
					cv::minMaxLoc(kernel, NULL, &max_val, NULL, &max_loc);
					if (max_loc.x == kernel_size && max_loc.y == kernel_size){
						peak_table.at<double>(j, i) = max_val;
						//cout << "find max = " << max_val << endl;
					}
				}
			}

		}
	}
	double first_peak = 0.0;
	double second_peak = 0.0;
	Point max_peak_loc;
	minMaxLoc(peak_table, NULL, &first_peak, NULL, &max_peak_loc);
	//cout << "loc.x = " << max_peak_loc.x << ",loc.y = " << max_peak_loc.y << endl;
	peak_table.at<double>(max_peak_loc.y, max_peak_loc.x) = 0.0;
	minMaxLoc(peak_table, NULL, &second_peak, NULL, &max_peak_loc);
	//cout << "loc.x = " << max_peak_loc.x << ",loc.y = " << max_peak_loc.y << endl;

	//�O��������O
	if (first_peak == 0.0){
		return(0.0);
	}
	if (first_peak != 0.0 && second_peak == 0.0){
		return(-100.0);
	}

	double sn = first_peak / second_peak;
	//cout << "sn = " << sn << endl;

	return(sn);
}

//�q�X�g�O�������Ƃ�A���10����臒l���擾����
int mw::makeHistgram(Mat src)
{
	Mat_<int> hist = Mat::zeros(1, 256, CV_32SC1);

	int tmpi;
	for (int j = 0; j < src.rows; j++){
		for (int i = 0; i < src.cols; i++){
			tmpi = (int)src.at<float>(j, i);
			if (0 <= tmpi && tmpi < 256){
				hist.at<int>(0, tmpi) = hist.at<int>(0, tmpi) + 1;
			}
		}
	}

	int sum = 0;
	int top10per = (int)(src.rows * src.cols * 0.1);
	int pixTop10per = 0;
	for (int cou = 255; cou > 0; cou--){
		sum += hist.at<int>(0, cou);

		if (top10per < sum){
			pixTop10per = cou;
			break;
		}
	}
	//makeGraph(hist);
	return(pixTop10per);
}

//void mw::plot(vector<double> x, vector<double> y, vector<vector<double> > z)
//{
//	CGnuplot gp;
//	gp.Plot(x, y, z);
//
//	__KEYWAIT__;
//}

//�n�������쐬
Mat_<float> mw::makeHannWindow(int width, int height)
{
	Mat_<float> window = Mat::zeros(height, width, CV_32FC1);
	
	float x = 0.0, y = 0.0;
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			x = (float)(i - width / 2.0) / width ;
			y = (float)(j - height / 2.0) / height;
			window.at<float>(j, i) = (float)((0.5 - 0.5* cos(2 * M_PI * x + M_PI))
				*( 0.5 - 0.5*cos(2 * M_PI * y + M_PI)));
		}
	}
	makeCsv<float>(window, "hann_window.csv");
	//dispMat<float>(window);
	//makeGraph3d(window);
	return(window);

}