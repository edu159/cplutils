import numpy as np
import itertools
import numpy as np
import re
import time
from cplutils.task import Task

class LogReader(Task):
    def __init__(self, fname, sections=[], columns=[{}]):
        Task.__init__(self)
        self.fname = fname
        self.sections = sections
        self.columns = columns
        self.output = {}

    #TODO: Add support for selecting certain columns by column number
    def _run(self):
        with open(self.fname, 'r') as log_file:
            lines = log_file.readlines()
        remaining_sections = [v for v in self.sections]
        cols = self.columns[0].copy()
        start = False
        run_steps = 0
        output_steps = 0
        initial_step = 0
        dump_every = 0
        s = 0
        sec_no = 0
        for sec in self.sections:
            self.output[sec] = {"time": [], "fields": {}}
            sec_no += 1
        current_section = None
        for ln, l in enumerate(lines):
            line_split = [v for v in l[:-1].split(" ") if v]
            # Avoid empty lines
            if not line_split:
                continue
            if current_section is None and self.sections:
                for sec in remaining_sections:
                    if "[%s]"%sec in l:
                        current_section = sec
                        remaining_sections.remove(current_section)
                        break
            else:
                if start:
                    for n, col in cols.items():
                        try:
                            # print line_split
                            self.output[current_section]["fields"][col][s] = float(line_split[1:][n])
                        except Exception as e:
                            raise Exception("Some fields are missing in line %d." % ln)
                    s += 1
                elif run_steps > 0 and line_split[0] == "Step":
                    start = True
                    try:
                        initial_step = int([v for v in lines[ln+1][:-1].split(" ") if v][0])
                        steps1 = int([v for v in lines[ln+2][:-1].split(" ") if v][0])
                        dump_every =  steps1 - initial_step
                        output_steps = run_steps / dump_every + 1
                        # By default if no sections provided. Return all the runs.
                        if not self.sections:
                            sec_name = "run%d" % sec_no
                            sec_no += 1
                            current_section = sec_name
                            self.output[sec_name] = {"time": [], "fields": {}}
                    except:
                        output_steps = 1

                    if not cols:
                        for n, col in enumerate(line_split[1:]):
                           cols[n] = col
                           self.output[current_section]["fields"][col] = np.zeros(output_steps)
                elif line_split[0] == "run":
                    try:
                        run_steps = int(line_split[1])
                    except:
                        pass
                # print s, output_steps
                if start and s == output_steps:
                    self.output[current_section]["time"] = xrange(output_steps)
                    start = False
                    current_section = None
                    run_steps = 0
                    output_steps = 0
                    initial_step = 0
                    dump_every = 0
                    s = 0
                    cols = {}
        # If the run has not been completed. End of log file.
        if current_section is not None and s < (output_steps):
            for c in self.output[current_section]["fields"].values():
                c.resize(s, refcheck=False)
                self.output[current_section]["time"] = xrange(s)

        return self.output


def log_reader(fname, sections=[], columns=[{}]):
    log = LogReader(fname, sections, columns)
    log.run()
    return log.output

class ChunkReader(Task):
    def __init__(self, fpath, chunktype, twrite=1.0, time_steps=[], cols=None, order="frame", zerostep=False):
        Task.__init__(self)
        self.fpath = fpath
        self.chunktype = chunktype
        self.twrite = twrite
        self.time_steps = time_steps
        self.cols = cols
        self.order = order
        self.output = {}
        self.zerostep = zerostep

    def _molchunk_frame_iter(self, arr_shape):
        return np.ndindex(arr_shape)

    def _molchunk_row_iter(self, arr_shape):
        index_iter = np.ndindex(arr_shape)
        while True:
            index = index_iter.next()
            yield (index[1], index[0])

    def _get_cells(self, frame0, dims=3):
        ncells = len(frame0)
        split_line = re.split("[ /\n]+", frame0[0].strip())
        xi0 = list([0.0]*3)
        xi1 = list([0.0]*3)
        xi0[0:dims] = map(float, split_line[1:dims+1])
        ncx = ncy = ncz = 1
        axis0 = [xi0[0]]
        axis1 = [xi0[1]]
        axis2 = [xi0[2]]
        # We assume no more than dims=3 (3D space)
        for l in frame0[1:]:
            split_line = re.split("[ /\n]+", l.strip())
            xi1[0:dims] = map(float, split_line[1:dims+1])
            if xi0[0] == xi1[0]:
                if xi0[1] == xi1[1]:
                    ncz += 1
                    axis2.append(xi1[2])
                else:
                    ncz = 1
                    axis2 = [axis2[0]]
                    ncy += 1
                    axis1.append(xi1[1])
                    xi0[1] = xi1[1]
            else:
                ncy = 1
                axis1 = [axis1[0]]
                ncz = 1
                axis2 = [axis2[0]]
                ncx += 1
                axis0.append(xi1[0])
                xi0[0] = xi1[0]
                xi0[1] = xi1[1]
        axis0 = np.array(axis0)
        axis1 = np.array(axis1)
        axis2 = np.array(axis2)
        return [ncx, ncy, ncz][0:dims], [axis0, axis1, axis2][0:dims]

    def _run(self):
        lines = []
        with open(self.fpath, 'r') as chunk_file:
            #TODO: Load up to certain number of frames depending on the max_RAM param
            lines = chunk_file.readlines()

        if self.chunktype == "molecule":
            skip_cols = 1
            header_lines = 3
            frame_header_lines = 1
            rows_per_frame = int(lines[3].split()[1])
            no_cols = len(lines[2].split()[skip_cols:-1])
        elif "spatial" in self.chunktype:
            dims = int(self.chunktype.split("-")[1][0])
            if dims not in [1, 2, 3]:
                raise Exception("Use spatial-1d, spatial-2d or spatial-3d.")
            skip_cols = dims + 1
            header_lines = 3
            frame_header_lines = 1
            rows_per_frame = int(lines[3].split()[1])
            no_cols = len(lines[2].split()[skip_cols:-1])
        elif self.chunktype == "trajectory-xyz":
            skip_cols = 1
            header_lines = 0
            frame_header_lines = 9
            rows_per_frame = int(lines[3].split()[0])
            no_cols = len(lines[8].split()[3:])
        elif self.chunktype == "ave/time-array":
            skip_cols = 1
            header_lines = 3
            frame_header_lines = 1
            rows_per_frame = int(lines[3].split()[1])
            no_cols = len(lines[2].split()[skip_cols:-1])

        liter = itertools.islice(lines, header_lines, None)
        no_lines = len(lines)
        if self.cols:
            no_cols = len(self.cols)
        if self.time_steps:
            no_frames = len(self.time_steps)
            time_steps_dt = (int(s) for s in np.array(self.time_steps)/(self.twrite))
        else:
            no_frames = (no_lines-header_lines) / (rows_per_frame + frame_header_lines)
            time_steps_dt = (s for s in xrange(0, no_frames))

        if "spatial" in self.chunktype:
            begin = header_lines + frame_header_lines
            end = begin + rows_per_frame
            ncells, axis = self._get_cells(lines[begin:end], dims=dims)
        zpad = 0
        if self.chunktype in ["molecule", "trajectory-xyz", "ave/time-array"]:
            if self.order == "frame": 
                if self.zerostep:
                    zpad = 1
                array_out2 = np.zeros((no_frames + zpad, rows_per_frame, no_cols), order='F')
                array_iter = self._molchunk_frame_iter(array_out2.shape[:-1])
            else:
                if self.zerostep:
                    zpad = 1
                array_out2 = np.zeros((rows_per_frame, no_frames + zpad, no_cols), order='C')
                array_iter = self._molchunk_row_iter((no_frames, rows_per_frame))
        elif "spatial" in self.chunktype:
            if self.zerostep:
                zpad = 1
            shape2 = tuple([no_frames + zpad] + ncells + [no_cols])
            shape = tuple([no_frames] + ncells + [no_cols])
            array_out2 = np.zeros(shape2, order='F')
            array_iter = self._molchunk_frame_iter(tuple(shape[:-1]))
        array_out = array_out2[zpad:,...]

        begin = 0
        end = 0
        frame =  time_steps_dt.next()
        iterator = enumerate(liter)
        step_index = {}
        if self.zerostep:
            step_index[0*self.twrite] = 0
            step_frame = 1
        else:
            step_frame = 0
        exit = False
        for n, l in iterator:
            if n == end:
                if exit: break
                begin =  frame*(frame_header_lines + rows_per_frame) + frame_header_lines
                end = begin + rows_per_frame
                iterator = itertools.islice(liter, begin, end)
                step_index[(frame+zpad)*self.twrite] = step_frame
                step_frame += 1
                try:
                    frame = time_steps_dt.next()
                except StopIteration:
                    exit = True

            else:
                if n >= begin and n < end:
                    split_line = np.fromstring(l, sep=' ')[skip_cols:]
                    selector = array_iter.next() + (slice(None),)
                    array_out[selector] = split_line


        self.output = {"data": array_out2, "tidx": step_index}
        if "spatial" in self.chunktype:
            self.output["axis"] = axis



def chunk_reader(fpath, chunktype, twrite=1.0, time_steps=[], cols=None, order="frame", zerostep=False):
    chunk = ChunkReader(fpath, chunktype, twrite, time_steps, cols, order, zerostep)
    chunk.run()
    return chunk.output
