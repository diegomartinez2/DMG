import cv2
import numpy as np

cap = cv2.VideoCapture(0)

while True:
    ret, frame = cap.read()
    if not ret:
        break

    h, w, _ = frame.shape
    noisy = frame.astype(np.float32)

    # 1. Ruido granular (estática)
    noise = np.random.normal(0, 28, (h, w, 1)).astype(np.float32)
    noisy = np.clip(noisy + noise, 0, 255).astype(np.uint8)

    # 2. Líneas horizontales de interferencia
    for _ in range(np.random.randint(2, 7)):
        y = np.random.randint(0, h - 1)
        band_h = np.random.randint(1, 5)
        shift = np.random.randint(-25, 25)
        noisy[y:y + band_h] = np.roll(noisy[y:y + band_h], shift, axis=1)

    cv2.imshow("Webcam con interferencia", noisy)

    if cv2.waitKey(1) & 0xFF == ord("q"):
        break

cap.release()
cv2.destroyAllWindows()
