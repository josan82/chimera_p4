import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 12, 41623])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwtOfYdVCWJhbGxTY2FsZXEDSwtHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLC0c/8AAAAAAAAH2HVQVjb2xvcnEFSwtLAH1xBihLAV1xB0sBYUsCXXEISwJhSwNdcQlLA2FLBF1xCksEYUsFXXELSwVhSwZdcQxLBmFLB11xDUsHYUsIXXEOSwhhSwldcQ9LCWFLCl1xEEsKYXWHVQpyaWJib25UeXBlcRFLC0sAfYdVCnN0aWNrU2NhbGVxEksLRz/wAAAAAAAAfYdVDG1tQ0lGSGVhZGVyc3ETXXEUKE5OTk5OTk5OTk5OZVUMYXJvbWF0aWNNb2RlcRVLC0sBfYdVCnZkd0RlbnNpdHlxFksLR0AUAAAAAAAAfYdVBmhpZGRlbnEXSwuJfYdVDWFyb21hdGljQ29sb3JxGEsLTn2HVQ9yaWJib25TbW9vdGhpbmdxGUsLSwB9h1UJYXV0b2NoYWlucRpLC4h9h1UKcGRiVmVyc2lvbnEbSwtLAH2HVQhvcHRpb25hbHEcfXEdVQhvcGVuZWRBc3EeiIlLCyhVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzFlM2dfbGlnYW5kLnNkZnEfVQtNREwgTU9ML1NERnEgTol0cSF9cSIoKFUwL2hvbWUvam9zZS9QaGFybWFjb3Bob3JlL0xpZ2FuZHMvMmFtYV9saWdhbmQuc2RmcSNoIE6JdHEkXXElSwRhKFU3L2hvbWUvam9zZS9QaGFybWFjb3Bob3JlL0xpZ2FuZHMvMmF4Nl9saWdhbmRfTmZpeGVkLnNkZnEmaCBOiXRxJ11xKEsGYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzF4cTNfbGlnYW5kLnNkZnEpaCBOiXRxKl1xK0sCYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzJhbWJfbGlnYW5kLnNkZnEsaCBOiXRxLV1xLksFYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzJhbTlfbGlnYW5kLnNkZnEvaCBOiXRxMF1xMUsDYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzFlM2tfbGlnYW5kLnNkZnEyaCBOiXRxM11xNEsBYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzRrN2FfbGlnYW5kLnNkZnE1aCBOiXRxNl1xN0sKYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzNybGpfbGlnYW5kLnNkZnE4aCBOiXRxOV1xOksIYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzNybGxfbGlnYW5kLnNkZnE7aCBOiXRxPF1xPUsJYShVMC9ob21lL2pvc2UvUGhhcm1hY29waG9yZS9MaWdhbmRzLzJodmNfbGlnYW5kLnNkZnE+aCBOiXRxP11xQEsHYXWHh3NVD2xvd2VyQ2FzZUNoYWluc3FBSwuJfYdVCWxpbmVXaWR0aHFCSwtHP/AAAAAAAAB9h1UPcmVzaWR1ZUxhYmVsUG9zcUNLC0sAfYdVBG5hbWVxREsLWA8AAAAySFZDX0xHRF9BXzIyMjZ9cUUoWA8AAAAxRTNLX1IxOF9BXzEwMDBdcUZLAWFYDwAAADFFM0dfUjE4X0FfMTAwMF1xR0sAYVgPAAAAMkFNQl8xN0hfQV8xMDAxXXFISwVhWAwAAAAzUkxMX1JMTF9BXzFdcUlLCWFYDAAAADNSTEpfUkxKX0FfMV1xSksIYVgPAAAAMkFNQV9ESFRfQV8xMDAxXXFLSwRhWA8AAAAxWFEzX1IxOF9BXzEwMDFdcUxLAmFYDwAAADJBTTlfVEVTX0FfMTAwMF1xTUsDYVgPAAAANEs3QV9ESFRfQV8xMDAxXXFOSwphWAwAAAAyQVg2X0hGVF9BXzFdcU9LBmF1h1UPYXJvbWF0aWNEaXNwbGF5cVBLC4l9h1UPcmliYm9uU3RpZmZuZXNzcVFLC0c/6ZmZmZmZmn2HVQpwZGJIZWFkZXJzcVJdcVMofXFUfXFVfXFWfXFXfXFYfXFZfXFafXFbfXFcfXFdfXFeZVUDaWRzcV9LC0sJSwCGfXFgKEsASwCGXXFhSwBhSwdLAIZdcWJLB2FLA0sAhl1xY0sDYUsISwCGXXFkSwhhSwZLAIZdcWVLBmFLAksAhl1xZksCYUsFSwCGXXFnSwVhSwpLAIZdcWhLCmFLAUsAhl1xaUsBYUsESwCGXXFqSwRhdYdVDnN1cmZhY2VPcGFjaXR5cWtLC0e/8AAAAAAAAH2HVRBhcm9tYXRpY0xpbmVUeXBlcWxLC0sCfYdVFHJpYmJvbkhpZGVzTWFpbmNoYWlucW1LC4h9h1UHZGlzcGxheXFuSwuIfYd1Lg=='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksLVQEgfYdVC2ZpbGxEaXNwbGF5cQNLC4l9h1UEbmFtZXEESwtYAwAAAFVOS32HVQVjaGFpbnEFSwtYAQAAACB9h1UOcmliYm9uRHJhd01vZGVxBksLSwJ9h1UCc3NxB0sLiYmGfYdVCG1vbGVjdWxlcQhLC0sAfXEJKEsBTl1xCksBSwGGcQthhksCTl1xDEsCSwGGcQ1hhksDTl1xDksDSwGGcQ9hhksETl1xEEsESwGGcRFhhksFTl1xEksFSwGGcRNhhksGTl1xFEsGSwGGcRVhhksHTl1xFksHSwGGcRdhhksITl1xGEsISwGGcRlhhksJTl1xGksJSwGGcRthhksKTl1xHEsKSwGGcR1hhnWHVQtyaWJib25Db2xvcnEeSwtOfYdVBWxhYmVscR9LC1gAAAAAfYdVCmxhYmVsQ29sb3JxIEsLTn2HVQhmaWxsTW9kZXEhSwtLAX2HVQVpc0hldHEiSwuJfYdVC2xhYmVsT2Zmc2V0cSNLC059h1UIcG9zaXRpb25xJF1xJShLAUsBhnEmSwFLAYZxJ0sBSwGGcShLAUsBhnEpSwFLAYZxKksBSwGGcStLAUsBhnEsSwFLAYZxLUsBSwGGcS5LAUsBhnEvSwFLAYZxMGVVDXJpYmJvbkRpc3BsYXlxMUsLiX2HVQhvcHRpb25hbHEyfVUEc3NJZHEzSwtK/////32HdS4='))
	atomInfo = cPickle.loads(base64.b64decode('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'))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECTRgBTn2HVQVhdG9tc3EDXXEEKF1xBShLFksXZV1xBihLFksfZV1xByhLF0sYZV1xCChLGEsZZV1xCShLGEspZV1xCihLGUsaZV1xCyhLGksbZV1xDChLGksfZV1xDShLG0scZV1xDihLHEsdZV1xDyhLHUseZV1xEChLHUsjZV1xEShLHksfZV1xEihLHksgZV1xEyhLIEshZV1xFChLIUsiZV1xFShLIksjZV1xFihLIksmZV1xFyhLIksnZV1xGChLI0skZV1xGShLJEslZV1xGihLJUsmZV1xGyhLJksoZV1xHChLJksqZV1xHShLK0ssZV1xHihLK0s0ZV1xHyhLLEstZV1xIChLLUsuZV1xIShLLUs+ZV1xIihLLksvZV1xIyhLL0swZV1xJChLL0s0ZV1xJShLMEsxZV1xJihLMUsyZV1xJyhLMkszZV1xKChLMks4ZV1xKShLM0s0ZV1xKihLM0s1ZV1xKyhLNUs2ZV1xLChLNks3ZV1xLShLN0s4ZV1xLihLN0s7ZV1xLyhLN0s8ZV1xMChLOEs5ZV1xMShLOUs6ZV1xMihLOks7ZV1xMyhLO0s9ZV1xNChLO0s/ZV1xNShLQEtBZV1xNihLQEtJZV1xNyhLQUtCZV1xOChLQktDZV1xOShLQktTZV1xOihLQ0tEZV1xOyhLREtFZV1xPChLREtJZV1xPShLRUtGZV1xPihLRktHZV1xPyhLR0tIZV1xQChLR0tNZV1xQShLSEtJZV1xQihLSEtKZV1xQyhLSktLZV1xRChLS0tMZV1xRShLTEtNZV1xRihLTEtQZV1xRyhLTEtRZV1xSChLTUtOZV1xSShLTktPZV1xSihLT0tQZV1xSyhLUEtSZV1xTChLUEtUZV1xTShLVUtWZV1xTihLVUtfZV1xTyhLVktXZV1xUChLV0tYZV1xUShLV0tZZV1xUihLWUtaZV1xUyhLWktbZV1xVChLWktfZV1xVShLW0tcZV1xVihLXEtdZV1xVyhLXUteZV1xWChLXUtjZV1xWShLXktfZV1xWihLXktgZV1xWyhLX0tpZV1xXChLYEthZV1xXShLYUtiZV1xXihLYktjZV1xXyhLYktmZV1xYChLYktoZV1xYShLY0tkZV1xYihLZEtlZV1xYyhLZUtmZV1xZChLZktnZV1xZShLaktrZV1xZihLakt0ZV1xZyhLa0tsZV1xaChLbEttZV1xaShLbEtuZV1xaihLbktvZV1xayhLb0twZV1xbChLb0t0ZV1xbShLcEtxZV1xbihLcUtyZV1xbyhLcktzZV1xcChLckt4ZV1xcShLc0t0ZV1xcihLc0t1ZV1xcyhLdEt+ZV1xdChLdUt2ZV1xdShLdkt3ZV1xdihLd0t4ZV1xdyhLd0t7ZV1xeChLd0t9ZV1xeShLeEt5ZV1xeihLeUt6ZV1xeyhLekt7ZV1xfChLe0t8ZV1xfShLf0uFZV1xfihLf0uGZV1xfyhLgEuGZV1xgChLgEuHZV1xgShLgUuCZV1xgihLgUuIZV1xgyhLgkuHZV1xhChLg0uEZV1xhShLg0uIZV1xhihLg0uJZV1xhyhLhEuKZV1xiChLhUuJZV1xiShLhkuUZV1xiihLh0uJZV1xiyhLiEuMZV1xjChLikuLZV1xjShLi0uMZV1xjihLi0uPZV1xjyhLi0uQZV1xkChLjEuNZV1xkShLjUuOZV1xkihLjkuPZV1xkyhLj0uTZV1xlChLj0uVZV1xlShLkEuRZV1xlihLkkuTZV1xlyhLlkuXZV1xmChLl0uYZV1xmShLl0uZZV1xmihLl0uaZV1xmyhLmkubZV1xnChLmkucZV1xnShLm0uiZV1xnihLnEudZV1xnyhLnEugZV1xoChLnUueZV1xoShLnUufZV1xoihLoEuhZV1xoyhLoUuiZV1xpChLokujZV1xpShLo0ukZV1xpihLpEulZV1xpyhLpEumZV1xqChLpkunZV1xqShLpkuoZV1xqihLpkupZV1xqyhLqkurZV1xrChLqku2ZV1xrShLq0vDZV1xrihLq0usZV1xryhLrEutZV1xsChLrUuxZV1xsShLrUu3ZV1xsihLrku3ZV1xsyhLr0u3ZV1xtChLsEu3ZV1xtShLsUu2ZV1xtihLsUuyZV1xtyhLskuzZV1xuChLs0u0ZV1xuShLs0u5ZV1xuihLtEu1ZV1xuyhLtUu2ZV1xvChLuEu5ZV1xvShLuEu+ZV1xvihLuUu6ZV1xvyhLukvCZV1xwChLu0u+ZV1xwShLvEu+ZV1xwihLvUu+ZV1xwyhLv0vCZV1xxChLwEvCZV1xxShLwUvCZV1xxihLxEvGZV1xxyhLxEvLZV1xyChLxUvNZV1xyShLxkvIZV1xyihLx0vNZV1xyyhLyEvKZV1xzChLyEvOZV1xzShLyUvNZV1xzihLykvMZV1xzyhLykvNZV1x0ChLy0vMZV1x0ShLy0vQZV1x0ihLzkvPZV1x0yhL0EvRZV1x1ChL0UvSZV1x1ShL0UvTZV1x1ihL00vUZV1x1yhL00vVZV1x2ChL00vWZV1x2ShL1kvXZV1x2ihL10vYZV1x2yhL2EvZZV1x3ChL2EvdZV1x3ShL2UvaZV1x3ihL2kvbZV1x3yhL20vcZV1x4ChL20veZV1x4ShL3EvdZV1x4ihL3kvfZV1x4yhL4Ev+ZV1x5ChL4UvoZV1x5ShL4kvpZV1x5ihL40v2ZV1x5yhL5Ev+ZV1x6ChL5Uv/ZV1x6ShL5kv/ZV1x6ihL50v/ZV1x6yhL6Ev4ZV1x7ChL6Uv5ZV1x7ShL6kvrZV1x7ihL6kvwZV1x7yhL60vxZV1x8ChL7EvtZV1x8ShL7Ev3ZV1x8ihL7Uv4ZV1x8yhL7kvvZV1x9ChL7kv5ZV1x9ShL70v6ZV1x9ihL8Ev8ZV1x9yhL8Uv9ZV1x+ChL8kv3ZV1x+ShL8kv7ZV1x+ihL80v1ZV1x+yhL80v+ZV1x/ChL9Ev2ZV1x/ShL9Ev3ZV1x/ihL9Uv6ZV1x/yhL9kv+ZV1yAAEAAChL+Ev7ZV1yAQEAAChL+Uv8ZV1yAgEAAChL+kv9ZV1yAwEAAChL+0v/ZV1yBAEAAChL/Ev9ZV1yBQEAAChNAAFNAQFlXXIGAQAAKE0AAU0KAWVdcgcBAAAoTQEBTQIBZV1yCAEAAChNAgFNAwFlXXIJAQAAKE0CAU0EAWVdcgoBAAAoTQQBTQUBZV1yCwEAAChNBQFNBgFlXXIMAQAAKE0FAU0KAWVdcg0BAAAoTQYBTQcBZV1yDgEAAChNBwFNCAFlXXIPAQAAKE0IAU0JAWVdchABAAAoTQgBTQ4BZV1yEQEAAChNCQFNCgFlXXISAQAAKE0JAU0LAWVdchMBAAAoTQoBTRQBZV1yFAEAAChNCwFNDAFlXXIVAQAAKE0MAU0NAWVdchYBAAAoTQ0BTQ4BZV1yFwEAAChNDQFNEQFlXXIYAQAAKE0NAU0TAWVdchkBAAAoTQ4BTQ8BZV1yGgEAAChNDwFNEAFlXXIbAQAAKE0QAU0RAWVdchwBAAAoTREBTRIBZWVVBWxhYmVsch0BAABNGAFYAAAAAH2HVQhoYWxmYm9uZHIeAQAATRgBiH2HVQZyYWRpdXNyHwEAAE0YAUc/yZmZoAAAAH2HVQtsYWJlbE9mZnNldHIgAQAATRgBTn2HVQhkcmF3TW9kZXIhAQAATRgBSwF9h1UIb3B0aW9uYWxyIgEAAH1VB2Rpc3BsYXlyIwEAAE0YAUsCfYd1Lg=='))
	crdInfo = cPickle.loads(base64.b64decode('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'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (7, 0, {}), 'display': (7, True, {}), 'id': (7, 100, {}), 'vrmlString': ['#VRML V2.0 utf8\nTransform {\n\ttranslation -3.083495 34.747514 4.504555\n\tchildren [\n\n\n\t\tShape {\n\t\t\tappearance Appearance {\n\t\t\t\tmaterial Material {\n\t\t\t\t\tambientIntensity 1\n\t\t\t\t\tdiffuseColor 0.000000 1.000000 1.000000\n\t\t\t\t\ttransparency 0.300000\n\t\t\t\t}\n\t\t\t}\n\t\t\tgeometry Sphere {\n\t\t\t\tradius 0.750000\n\t\t\t}\n\t\t}\n\t]\n}', '#VRML V2.0 utf8\nTransform {\n\ttranslation -3.682074 35.054971 4.835725\n\tchildren [\n\t\tTransform {\n\t\t\trotation 0 0 1 -2.045300\n\t\t\tchildren [\n\t\t\t\tTransform {\n\t\t\t\t\trotation 1 0 0 -0.457336\n\t\t\t\t\tchildren [\n\n\n\t\t\t\t\t\tShape {\n\t\t\t\t\t\t\tappearance Appearance {\n\t\t\t\t\t\t\t\tmaterial Material {\n\t\t\t\t\t\t\t\t\tambientIntensity 1\n\t\t\t\t\t\t\t\t\tdiffuseColor 0.000000 1.000000 1.000000\n\t\t\t\t\t\t\t\t\ttransparency 0.500000\n\t\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\tgeometry Cone {\n\t\t\t\t\t\t\t\tbottomRadius 0.750000\n\t\t\t\t\t\t\t\theight 1.500000\n\t\t\t\t\t\t\t\tbottom FALSE\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t}\n\t\t\t\t\t]\n\t\t\t\t}\n\t\t\t]\n\t\t}\n\t]\n}',
'#VRML V2.0 utf8\nTransform {\n\ttranslation 3.197117 26.753558 4.342462\n\tchildren [\n\n\n\t\tShape {\n\t\t\tappearance Appearance {\n\t\t\t\tmaterial Material {\n\t\t\t\t\tambientIntensity 1\n\t\t\t\t\tdiffuseColor 1.000000 0.000000 1.000000\n\t\t\t\t\ttransparency 0.300000\n\t\t\t\t}\n\t\t\t}\n\t\t\tgeometry Sphere {\n\t\t\t\tradius 0.750000\n\t\t\t}\n\t\t}\n\t]\n}', '#VRML V2.0 utf8\nTransform {\n\ttranslation -3.083495 34.747514 4.504555\n\tchildren [\n\n\n\t\tShape {\n\t\t\tappearance Appearance {\n\t\t\t\tmaterial Material {\n\t\t\t\t\tambientIntensity 1\n\t\t\t\t\tdiffuseColor 1.000000 0.000000 1.000000\n\t\t\t\t\ttransparency 0.300000\n\t\t\t\t}\n\t\t\t}\n\t\t\tgeometry Sphere {\n\t\t\t\tradius 0.750000\n\t\t\t}\n\t\t}\n\t]\n}', '#VRML V2.0 utf8\nTransform {\n\ttranslation -3.682074 35.054971 4.835725\n\tchildren [\n\t\tTransform {\n\t\t\trotation 0 0 1 -2.045300\n\t\t\tchildren [\n\t\t\t\tTransform {\n\t\t\t\t\trotation 1 0 0 -0.457336\n\t\t\t\t\tchildren [\n\n\n\t\t\t\t\t\tShape {\n\t\t\t\t\t\t\tappearance Appearance {\n\t\t\t\t\t\t\t\tmaterial Material {\n\t\t\t\t\t\t\t\t\tambientIntensity 1\n\t\t\t\t\t\t\t\t\tdiffuseColor 1.000000 0.000000 1.000000\n\t\t\t\t\t\t\t\t\ttransparency 0.500000\n\t\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\tgeometry Cone {\n\t\t\t\t\t\t\t\tbottomRadius 0.750000\n\t\t\t\t\t\t\t\theight 1.500000\n\t\t\t\t\t\t\t\tbottom FALSE\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t}\n\t\t\t\t\t]\n\t\t\t\t}\n\t\t\t]\n\t\t}\n\t]\n}',
'#VRML V2.0 utf8\nTransform {\n\ttranslation 1.548304 30.725067 4.859316\n\tchildren [\n\n\n\t\tShape {\n\t\t\tappearance Appearance {\n\t\t\t\tmaterial Material {\n\t\t\t\t\tambientIntensity 1\n\t\t\t\t\tdiffuseColor 0.500000 0.250000 0.000000\n\t\t\t\t\ttransparency 0.300000\n\t\t\t\t}\n\t\t\t}\n\t\t\tgeometry Sphere {\n\t\t\t\tradius 0.750000\n\t\t\t}\n\t\t}\n\t]\n}', '#VRML V2.0 utf8\nTransform {\n\ttranslation -0.654114 31.771784 4.440472\n\tchildren [\n\n\n\t\tShape {\n\t\t\tappearance Appearance {\n\t\t\t\tmaterial Material {\n\t\t\t\t\tambientIntensity 1\n\t\t\t\t\tdiffuseColor 0.500000 0.250000 0.000000\n\t\t\t\t\ttransparency 0.300000\n\t\t\t\t}\n\t\t\t}\n\t\t\tgeometry Sphere {\n\t\t\t\tradius 0.750000\n\t\t\t}\n\t\t}\n\t]\n}'], 'name': (7, u'p4', {})}
	colors = {u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'gold': ((1, 0.843137, 0), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'),
u'Pr': ((0.85098, 1, 0.780392), 1, u'default'), u'deep pink': ((1, 0.0784314, 0.576471), 1, u'default'), u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), u'Pu': ((0, 0.419608, 1), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'purple': ((0.627451, 0.12549, 0.941176), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'),
u'Mt': ((0.921569, 0, 0.14902), 1, u'default'), u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), u'sky blue': ((0.529412, 0.807843, 0.921569), 1, u'default'), u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'plum': ((0.866667, 0.627451, 0.866667), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'dodger blue': ((0.117647, 0.564706, 1), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'),
u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'light green': ((0.564706, 0.933333, 0.564706), 1, u'default'),
u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), u'Tm': ((0, 0.831373, 0.321569), 1, u'default'), u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), u'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'salmon': ((0.980392, 0.501961, 0.447059), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'),
u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'), u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), u'light gray': ((0.827451, 0.827451, 0.827451), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 14, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (16, (u'deep pink', (1, 0.0784314, 0.576471, 1)), {(u'green', (0, 1, 0, 1)): [15], (u'light green', (0.564706, 0.933333, 0.564706, 1)): [3], (u'dodger blue', (0.117647, 0.564706, 1, 1)): [8], (u'', (0.358033, 0.260402, 0.804281, 1)): [10], (u'N', (0.188235, 0.313725, 0.972549, 1)): [13], (u'F', (0.564706, 0.878431, 0.313725, 1)): [12], (u'purple', (0.627451, 0.12549, 0.941176, 1)): [9], (u'gold', (1, 0.843137, 0, 1)): [7], (u'sky blue', (0.529412, 0.807843, 0.921569, 1)): [1], (u'tan', (0.823529, 0.705882, 0.54902, 1)): [0], (u'O', (1, 0.0509804, 0.0509804, 1)): [11], (u'plum', (0.866667, 0.627451, 0.866667, 1)): [2], (u'light gray', (0.827451, 0.827451, 0.827451, 1)): [5], (u'salmon', (0.980392, 0.501961, 0.447059, 1)): [4], (u'yellow', (1, 1, 0, 1)): [14]})
	viewerInfo = {'cameraAttrs': {'center': (0.10849997615814, 30.985500019073, 4.3005), 'fieldOfView': 25.350031531442, 'nearFar': (15.370701651827, -9.9165740527803), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 4.3005}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 10.031640289183, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 1, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 15, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(11, 'Chimera default', 'rounded', u'unknown'), (12, 'Chimera default', 'rounded', u'unknown'), (13, 'Chimera default', 'rounded', u'unknown'), (14, 'Chimera default', 'rounded', u'unknown'), (15, 'Chimera default', 'rounded', u'unknown'), (16, 'Chimera default', 'rounded', u'unknown'), (17, 'Chimera default', 'rounded', u'unknown'), (18, 'Chimera default', 'rounded', u'unknown'), (19, 'Chimera default', 'rounded', u'unknown'), (20, 'Chimera default', 'rounded', u'unknown'), (21, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.35740674433659325, 0.6604015517481454, -0.6604015517481455), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.35740674433659325, 0.6604015517481454, 0.6604015517481455), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.25056280708573153, 0.25056280708573153, 0.9351131265310293), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (721, 563)
	xformMap = {0: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 1: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 2: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 3: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 4: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 5: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 6: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 7: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True),
8: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 9: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 10: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 557: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 558: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 559: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 560: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 561: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True),
562: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True), 563: (((-0.42407096693773, 0.81710586963004, -0.3905148048623), 24.980823802706), (-4.4111083822129, 0.53368823181186, 6.7856007597191), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 6: True, 7: True, 8: True, 9: True, 10: True, 557: True, 558: True, 559: True, 560: True, 561: True, 562: True, 563: True, 564: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

