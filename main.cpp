#include <Novice.h>
#include "GeometryUtility.h"
#include <corecrt_math_defines.h>

const char kWindowTitle[] = "学籍番号";

static const int kRowHeight = 20;
static const int kColumnWidth = 60;

GeometryUtility geometryUtility;

Matrix4x4 DirectionToDirection(const Vector3& from, const Vector3& to) {
	
	Vector3 fromNormalized = geometryUtility.normalize(from);
	Vector3 toNormalized = geometryUtility.normalize(to);

	
	Vector3 axis = geometryUtility.cross(fromNormalized, toNormalized);
	float dot = geometryUtility.Dot(fromNormalized, toNormalized);

	
	if (dot > 0.9999f) {
		
		return {
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		};
	} else if (dot < -0.9999f) {
		
		Vector3 orthoAxis = geometryUtility.Perpendicular(fromNormalized);
		return geometryUtility.MakeRotateAxisAngle(orthoAxis, float(M_PI));
	}

	axis = geometryUtility.normalize(axis);

	float angle = acosf(dot);

	return geometryUtility.MakeRotateAxisAngle(axis, angle);

}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 500, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		Vector3 from0 = geometryUtility.normalize(Vector3{ 1.0f,0.7f,0.5f });
		Vector3 to0 = geometryUtility.Multiply(-1.0f, from0);
		Vector3 from1 = geometryUtility.normalize(Vector3{ -0.6f,0.9f,0.2f });
		Vector3 to1 = geometryUtility.normalize(Vector3{ 0.4f,0.7f,-0.5f });
		


		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		Matrix4x4 rotateMatrix0 = DirectionToDirection(
			geometryUtility.normalize(Vector3{ 1.0f,0.0f,0.0f }), geometryUtility.normalize(Vector3{ -1.0f,0.0f,0.0f }));
		Matrix4x4 rotateMatrix1 = DirectionToDirection(from0, to0);
		Matrix4x4 rotateMatrix2 = DirectionToDirection(from1, to1);
		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		geometryUtility.MatrixScreenPrintf(0, 0, rotateMatrix0, "rotateMatrix0");
		geometryUtility.MatrixScreenPrintf(0, kRowHeight * 5 , rotateMatrix1, "rotateMatrix1");
		geometryUtility.MatrixScreenPrintf(0, kRowHeight * 10, rotateMatrix2, "rotateMatrix2");

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
