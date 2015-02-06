#pragma once

namespace mw
{
	void String2string(System::String ^ s, std::string& os);
	System::String^ String2string(std::string str);
	System::String^ double2String(double tmpd);
	System::String^ int2String(int tmpi);
	System::String^ getFilePath();
	System::String^ getFileListPath();
	array<System::String^>^ getFileList();
	void MarshalString(System::String ^ s, string& os);
	string int2string(int number);

	Mat& csv2mat(string file_name);
	Mat_<float> csv2matF(string file_name);

	void makeGraph2d(cv::Mat hist);
	void makeGraph3d(cv::Mat src);

	void swqp(int& x, int& y);
	void makeRandom(int* ans);

	Point calGlavityCenter(Mat_<double> src);
	double calZNCC(Mat_<int>& bf, Mat_<int>& af);
	double calNMIinMask(Mat_<int> bf, Mat_<int> af, Mat_<int> mask, int BIN_DOWN);
	double calSN(Mat_<double> src);
	double calSN_ByPeak(Mat_<double>& src);
	double calSN_ByPeak(Mat_<double>& src, int mountain_r);
	int makeHistgram(Mat src);
	//void dispMat(Mat_<int> src);
	Mat_<float> makeHannWindow(int width, int height);

	template<typename T_n>
	double calAve(Mat_<T_n> src)
	{
		int num_all_pixel = src.rows * src.cols;
		T_n sum = 0.0;
		for (int j = 0; j < src.rows; j++){
			for (int i = 0; i < src.cols; i++){
				sum += src.at<T_n>(j, i);
			}
		}
		double dst;
		dst = (double)sum / (double)num_all_pixel;
		return(dst);
	}

	////csvファイルを作るよ
	//T_nはsrcの型
	template<typename T_n>
	void makeCsv(cv::Mat src, char file_name[])
	{
		FILE *fp;

		errno_t error;
		if (error = fopen_s(&fp, file_name, "w") != 0){
			fprintf(stderr, "ファイルのオープンに失敗しました．\n");
		}

		double tmp;
		for (int j = 0; j < src.rows; ++j){
			for (int i = 0; i < src.cols; ++i){
				tmp = (double)src.at<T_n>(j, i);
				fprintf(fp, "%lf", tmp);
				if (i == src.cols - 1){
					fprintf(fp, "\n", tmp);
				}
				else{
					fprintf(fp, ",");
				}
			}
		}
		fclose(fp);
		cout << file_name << "のファイルを作成しました！";
	}




	template<class T_in , class T_out>
	Mat_<T_out> array_2_mat(array<T_in, 2>^ src)
	{
		Mat_<T_out> tmp_mat = cv::Mat_<T_out>::zeros(src->GetLength(0), src->GetLength(1));
		for (int j = 0; j < src->GetLength(0); j++){
			for (int i = 0; i < src->GetLength(1); i++){
				tmp_mat.at<T_out>(j, i) = (T_out)src[j, i];
			}
		}
		return(tmp_mat);
	}

	template<class T_in, class T_out>
	array<T_out, 2>^ mat_2_array(Mat_<T_in> src)
	{
		array<T_out, 2>^ dst = gcnew array<int, 2>(src.rows , src.cols);
		for (int j = 0; j < src.rows; j++){
			for (int i = 0; i < src.cols; i++){
				dst[j, i] = (T_out)src.at<T_in>(j, i);
			}
		}
		return(dst);
	}

	template<class T_in>
	Mat array_2_ViewMat(array<T_in, 2>^ src)
	{
		Mat dst = cv::Mat(src->GetLength(0), src->GetLength(1), CV_8UC3, Scalar(0, 0, 0));
		int tmpi;
		for (int j = 0; j < src->GetLength(0); j++){
			for (int i = 0; i < src->GetLength(1); i++){
				tmpi = (uchar)src[j, i];
				if (tmpi < 0){ tmpi = 0; };
				if (tmpi > 254){ tmpi = 254; };

				dst.at<Vec3b>(j,i) = Vec3b(tmpi, tmpi, tmpi);
			}
		}
		return(dst);
	}

	//Matの確認用
	template<class T_n>
	void dispMat(Mat src)
	{
		Mat_<double> pic_result = Mat_<double>::zeros(src.rows, src.cols);

		double max_val, min_val;
		cv::minMaxLoc(src, &min_val, &max_val, NULL, NULL);
		for (int j = 0; j < src.rows; j++){
			for (int i = 0; i < src.cols; i++){
				pic_result[j][i] = (src.at<T_n>(Point(i, j)) - min_val) / (max_val - min_val);
			}
		}
		cv::namedWindow("一時的確認", CV_WINDOW_AUTOSIZE | CV_WINDOW_FREERATIO);
		cv::imshow("一時的確認", pic_result);
		waitKey(0);
	}

	//Matの確認用
	template<class T_n>
	void dispMat(Mat src , int number)
	{
		Mat_<double> pic_result = Mat_<double>::zeros(src.rows, src.cols);

		double max_val, min_val;
		cv::minMaxLoc(src, &min_val, &max_val, NULL, NULL);
		for (int j = 0; j < src.rows; j++){
			for (int i = 0; i < src.cols; i++){
				pic_result[j][i] = (src.at<T_n>(Point(i, j)) - min_val) / (max_val - min_val);
			}
		}
		cv::namedWindow(mw::int2string(number), CV_WINDOW_AUTOSIZE | CV_WINDOW_FREERATIO);
		cv::imshow(mw::int2string(number), pic_result);
		//waitKey(0);
	}

};
