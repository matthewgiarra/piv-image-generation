function img = getSensorNoise(Camera)
	
    img = abs(Camera.SensorNoiseStd * randn(Camera.PixelRows, Camera.PixelColumns));

end