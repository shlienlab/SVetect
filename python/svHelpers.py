"""
Helper functions commonly used in svetect modules.

"""

import logging
import warnings
import yaml

def configure(obj, config, class_name='detexy'):

    """
    python3 configure(obj, config, class_name)

    Arguments:

        - obj: object to configure
        - config: path to config file
        - class_name: section of config containing class-specific properties

    Overwrites any default properties with those passed in the config.

    """

    with open(config, 'r') as f:
        options = yaml.safe_load(f)
        options = options.get(class_name, {})
        PYSAM_ZERO_BASED = options.get('PYSAM_ZERO_BASED', 1)

    obj.PYSAM_ZERO_BASED = PYSAM_ZERO_BASED
    for option in options:
        if isinstance(options[option], dict):
            default = getattr(obj, option)
            nested = options[option]
            for opt in nested.keys():
                default_arg = default.get(opt, 'MISSING')
                if default_arg is 'MISSING':
                    warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                    continue
                if not isinstance(nested[opt], type(default_arg)) and default_arg is not None:
                    warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(nested[opt])))
                    continue
                default[opt] = nested[opt]
            setattr(obj, option, default)
        else:
            default_arg = getattr(obj, option, 'MISSING')
            if default_arg is 'MISSING':
                warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                continue
            if not isinstance(options[option], type(default_arg)) and default_arg is not None:
                warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(options[option])))
                continue
            setattr(obj, option, options[option])

def setup_log(log_id, logfile, level=logging.INFO, mode='w'):

    """
    Set up log file.

    """

    logging.basicConfig(level=level)
    logger = logging.getLogger(log_id)
    logger.setLevel(level)
    handler = logging.FileHandler(logfile, mode=mode)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def clean_chromosomes(chromosomes):

    """
    Ensure that chromosome is of type str.

    Valid chromosomes are '0' to '22' (inclusive) and 'X' and 'Y'.

    """

    return [str(int(x)) if x in list(range(0, 23)) else str(x) for x in chromosomes]

