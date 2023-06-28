import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from manas_cafa5.model import Model
from manas_cafa5.protein import Protein
from manas_cafa5.protein_cache import ProteinCache

file = 'data/go_terms_train_set_maxlen500_minmembers75.tsv'

m = Model.cnn1d([],1024)
trainingset = Model.parse_trainingset(file)
graph = Protein.build_graph('data/go-basic.obo')

protein_cache = ProteinCache('data')
loaded = m.compile(trainingset, graph, 'P68510', protein_cache)

checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
    filepath=os.path.join('ch', "ckpt_{epoch}"),
    save_weights_only=True,
)
loaded.fit(epochs=5, callbacks=[checkpoint_callback])
