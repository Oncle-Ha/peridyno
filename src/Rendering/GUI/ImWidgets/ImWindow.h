#pragma once

#include <vector>
#include <memory>

namespace dyno
{
	struct Picture;
	class RenderEngine;
	class SceneGraph;

	class ImWindow
	{
	public:
		void initialize(float scale);
		void draw(RenderEngine* engine, SceneGraph* scene);

		bool cameraLocked();
		bool saveScreen();
		int screenTime() { return mScreenTime;}
	private:
		bool mDisenableCamera = false;
		bool mSaveScreen = false;
		int mScreenTime = 0;
		std::vector<std::shared_ptr<Picture>> mPics;
	};
}
