function img = getSensorNoise(Camera)
	
    img = abs(Camera.SensorNoiseMean + Camera.SensorNoiseStd * randn(Camera.PixelRows, Camera.PixelColumns));

end