# import analysisUtils as au
"""
# from astroquery.ipac.ned import Ned
# result_table = Ned.query_object("MCG12-02-001")
# z = result_table['Redshift'].data.data
"""
#redshift for some sources.
z_d = {'VV705': 0.04019,
       'UGC5101':0.03937,
       'UGC8696':0.03734,
       'MCG12' : 0.015698,
       'VV250':0.03106}


# def print_logger_header(title, logger):
#     separator = "-" * len(title)
#     logger.info(separator)
#     logger.info(title)
#     logger.info(separator)

def print_logger_header(title, logger):
    width = len(title) + 4  # Add padding to the width
    top_bottom_border = "+" + "-" * width + "+"
    side_border = "| " + title + " |"

    logger.info(top_bottom_border)
    logger.info(side_border)
    logger.info(top_bottom_border)



def deprecated(old_name, new_name):
    def decorator(func):
        def wrapper(*args, **kwargs):
            warnings.warn(f"'{old_name}' is deprecated and "
                          f"will be removed in a future version. "
                          f"Use '{new_name}' instead.",
                          category=DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)
        return wrapper
    return decorator


