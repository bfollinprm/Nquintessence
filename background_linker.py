import background

def get_background_hubble(*args, **kwargs):
	background_class_instance = background.background(args, kwargs)
	hubble = background_class_instance.hubble

	return hubble