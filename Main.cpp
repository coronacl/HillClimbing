#include "stdafx.h"

//一次元ベクトル＜double>の最大値と最小値を取得
void calMax(vector<double> src , double& _max_val , int& _max_pos)
{
	if (src.size() == 0){
		cout << "vectorの中は空っぽだよ" << endl;
		return;
	}

	double max_val = src[0];
	int max_pos = 0;
	for (int depth = 1; depth < (int)src.size(); depth++){
		if (src[depth] > max_val){
			max_val = src[depth];
			max_pos = depth;
		}
	}
	_max_val = max_val;
	_max_pos = max_pos;
	//cout << "max_val = " << max_val << endl;
	//cout << "max_pos = " << max_pos << endl;
}

////山登り法
//data	探索対象
//pos_start	探索開始地点
//goal	探索終了地点
//glal_val	探索終了地点の大きさ
void hillClimbing(Mat_<double> data , Point pos_start , Point& goal , double& goal_val)
{
	Point search = pos_start;
	cout << "start_point =" << pos_start << endl;

	Point direction[5];
	direction[0] = Point(0, 0);
	direction[1] = Point(1, 0);
	direction[2] = Point(0, 1);
	direction[3] = Point(-1, 0);
	direction[4] = Point(0, -1);

	vector<double> altitude;
	altitude.resize(9);

	int alpha = 4;
	bool flag_search = true;
	int counter = 0;
	int tmpX, tmpY;
	double max_val;
	while (flag_search){
		if (counter++ > 100){
			break;
		}

		//各方向の高さを取得
		for (int i = 0; i < 5; i++){
			altitude[i] = 0.0;
			tmpX = search.x + alpha * direction[i].x;
			tmpY = search.y + alpha * direction[i].y;
			if (0 <= tmpX || tmpX < data.cols || 0 <= tmpY || tmpY < data.rows){
				altitude[i] = data.at<double>(tmpY, tmpX);
			}
		}

		////最も勾配が高い方向へ移動
		int max_index = 0;
		calMax(altitude, max_val, max_index);
		if (max_index == 0){
			if (alpha == 1){
				flag_search = false;
				break;
			}
			alpha = 1;
		}
		search.x = search.x + alpha * direction[max_index].x;
		search.y = search.y + alpha * direction[max_index].y;
	}
	cout << "  counter:" << counter << endl;
	cout << "  "<< search.x << "," << search.y << "::" << max_val << endl;

	goal = search;
	goal_val = max_val;
}

void calMaxByHillClimbing(Mat_<double> data)
{
	//初期点を作成
	vector<Point> list;
	for (int j = 1; j < 4; j++){
		for (int i = 1; i < 4; i++){
			list.push_back(Point((int)(data.cols / 4 * i), (int)(data.rows / 4 * j)));
		}
	}

	//各初期点を山登り法で探索
	vector<double> list_altitude;
	list_altitude.resize(list.size());
	vector<Point> list_goal_pos;
	list_goal_pos.resize(list.size());
	for (int i = 0; i < (int)list.size(); i++){
		hillClimbing(data, list[i] , list_goal_pos[i] , list_altitude[i]);
	}

	//最も高い山を返却
	int index;
	double tmpd;
	calMax(list_altitude,tmpd,index);
	cout << "最も高い山の位置は、" << list_goal_pos[index] << endl;
}

int main(array<System::String ^> ^args)
{
	Mat_<double> mountain = mw::csv2matF("mountain.csv");
	mw::makeGraph3d(mountain);

	double max_value;
	Point max_location;
	minMaxLoc(mountain, NULL, &max_value, NULL, &max_location);
	cout << "max_value = " << max_value << endl;
	cout << "max_location = " << max_location << endl << endl;
	
	calMaxByHillClimbing(mountain);

	int a; cin >> (int)a;//グラフを表示させ続けるため
	return(0);
}
