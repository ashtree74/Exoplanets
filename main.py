# The application of Kalman Filtering to the analysis
# of baseline 1D star photometry data (exoplanets transit)
# by Adam Jesionkiewicz (adam@jesion.pl)

import plotly.plotly
from plotly.graph_objs import *
import datetime, math
import logging

# Filter parameters (def: 0.008, 0.1)
FILTER_R = 0.002
FILTER_Q = 0.1
TIME_CRT = 0.000

FILENAME = 'data/data.txt'


class DataStream():
    """
    A few method to operate on raw data set.
    Purpose: convert data from "2458126.22243 -0.09398 0.0029" to:
    {'dif': 0.0029, 'dt': datetime.datetime(2018, 1, 7, 17, 20, 17, 952007),
    'ts': 2458126.22243, 'mag': -0.09398, 'f_mag': -0.09398}
    """
    cov = float('nan')
    x = float('nan')

    def __init__(self, file_name):
        """
        Constructor with automatic data loader
        """
        self.file_name = file_name

        self.A = 1
        self.B = 0
        self.C = 1
        self.R = FILTER_R
        self.Q = FILTER_Q

        self.data_stream = []

    def load_data(self):
        """
        Open and parse a local disk file (raw data)
        :return: a list of raw records
        """
        logging.debug('Loading data from file ({})...'.format(self.file_name))
        parsed_data = list()
        with open(self.file_name) as file_data:
            for line in file_data.readlines():
                temp = dict()
                if 'JD' in line:
                    continue
                line = line.split()
                temp['ts'], temp['mag'], temp['dif'] = float(line[0][:14]), float(line[1]), float(line[2])
                temp['f_mag'] = self.kalman_filter(temp['mag'])
                temp['dt'] = self.jd_to_datetime(temp['ts'])
                temp['dt_cor'] = self.jd_to_datetime(temp['ts'] - TIME_CRT)
                parsed_data.append(temp)
        logging.debug('  {} records loaded.'.format(len(parsed_data)))
        logging.debug(parsed_data[0])
        self.data_stream = parsed_data

    def kalman_filter(self, measurement):
        """
        Filters a measurement. More: https://en.wikipedia.org/wiki/Kalman_filter
        :param measurement: The measurement value to be filtered
        :return: The filtered value
        """
        u = 0
        if math.isnan(self.x):
            self.x = (1 / self.C) * measurement
            self.cov = (1 / self.C) * self.Q * (1 / self.C)
        else:
            pred_x = (self.A * self.x) + (self.B * u)
            pred_cov = ((self.A * self.cov) * self.A) + self.R

            # Kalman Gain
            k = pred_cov * self.C * (1 / ((self.C * pred_cov * self.C) + self.Q))

            # Correction
            self.x = pred_x + k * (measurement - (self.C * pred_x))
            self.cov = pred_cov - (k * self.C * pred_cov)

        return self.x

    def jd_to_datetime(self, jd):
        """
        Convert a Julian Day to a gregorian calendar.
        :param jd: example: 2458126.24171
        :return: datetime object
        """
        jd = jd + 0.5
        F, I = math.modf(jd)
        I = int(I)
        A = math.trunc((I - 1867216.25) / 36524.25)
        if I > 2299160:
            B = I + 1 + A - math.trunc(A / 4.)
        else:
            B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
        if G < 13.5:
            month = G - 1
        else:
            month = G - 13
        if month > 2.5:
            year = D - 4716
        else:
            year = D - 4715

        frac_days, day = math.modf(day)
        day = int(day)

        hours = frac_days * 24.
        hours, hour = math.modf(hours)

        mins = hours * 60.
        mins, min = math.modf(mins)

        secs = mins * 60.
        secs, sec = math.modf(secs)

        micro = round(secs * 1.e6)

        return datetime.datetime(year, month, day, int(hour), int(min), int(sec), int(micro))

    def get_data(self):
        if self.data_stream:
            return self.data_stream
        else:
            self.load_data()
            return self.data_stream


class VisualizeData():
    """
    Plot data on the diagram (plot.ly)
    """
    def __init__(self, data):
        self.data = data
        logging.debug('VisualizedData - initialize')

    def plot_graph(self, dataset):
        """
        Plot diagram on the server (plotly)
        :param dataset: choose 'mag', 'f_mag', or 'dif'
        """
        data = self.data
        diagrams = []

        for time_stamp, data_tag in dataset:
            data_x, data_y = [], []
            for item in data:
                data_x.append(item[time_stamp])
                data_y.append(item[data_tag])
            diagrams.append(Scatter(x=data_x, y=data_y, mode='markers'))

        layout = plotly.graph_objs.Layout(yaxis=dict(autorange='reversed'))
        data = Data(diagrams)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        plotly.plotly.plot(fig, filename='exo-line')


if __name__ == '__main__':

    logging.basicConfig(
        filename="debug.log",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s:%(message)s"
    )

    logging.debug('File name: {}'.format(FILENAME))

    app = DataStream(FILENAME)
    VisualizeData(app.get_data()).plot_graph([('dt','mag'), ('dt_cor','f_mag')])

